#include <algorithm>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "gdal.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "cpl_error.h"
#include "ogr_feature.h"
#include "ogr_geometry.h"
#include "ogr_core.h"
#include "ogr_spatialref.h"
#include "cpl_port.h"
#include "cpl_string.h"

using GeoTransform = std::array<double, 6>;
using ArrayIndex_t = uint32_t;
using Row = ArrayIndex_t;
using Column = ArrayIndex_t;
using ArrayCoordinate = std::pair<Column, Row>;
using Latitude = double;
using Longitude = double;
using GeographicCoordinate = std::pair<Longitude, Latitude>;

struct Arguments {
public:
    size_t raster_band;
    size_t cache_size;
    std::string point_layer;
    std::string input_raster;
    std::string input_points;
    std::string output_points;
    std::string output_format;

    Arguments(): raster_band(1),
                 cache_size(64),
                 point_layer(""),
                 input_raster(""),
                 input_points(""),
                 output_points(""),
                 output_format("GPKG") {};
    Arguments(const size_t raster_band,
              const size_t cache_size,
              const std::string point_layer,
              const std::string input_raster,
              const std::string input_points,
              const std::string output_points,
              const std::string output_format): raster_band(raster_band),
                                                cache_size(cache_size),
                                                point_layer(point_layer),
                                                input_raster(input_raster),
                                                input_points(input_points),
                                                output_points(output_points),
                                                output_format(output_format) {};
    ~Arguments() = default;
};

struct RasterShape {
public:
    size_t blockxsize;
    size_t blockysize;
    size_t size_x;
    size_t size_y;
    size_t num_blocks_x;
    size_t num_blocks_y;

    RasterShape() = default;
    ~RasterShape() = default;
};

ArrayCoordinate projection_to_rowcol(const GeoTransform& gt, const GeographicCoordinate& point) {
    return std::make_pair(
        static_cast<ArrayIndex_t>((point.first - gt[0]) / gt[1]),
        static_cast<ArrayIndex_t>((point.second - gt[3]) / gt[5])
    );
}

std::vector<ArrayCoordinate> projection_to_rowcol(const GeoTransform& gt, const std::vector<GeographicCoordinate>& points) {
    std::vector<ArrayCoordinate> output{};
    output.reserve(points.size());
    
    std::transform(points.begin(), points.end(), std::back_inserter(output), [&gt](const auto& elem) {
        return projection_to_rowcol(gt, elem);
    });

    return output;
}

GeographicCoordinate rowcol_to_projection(const GeoTransform& gt, const ArrayCoordinate& point) {
    return std::make_pair(
        gt[1] * point.first + gt[2] * point.second + gt[1] * 0.5 + gt[2] * 0.5 + gt[0],
        gt[4] * point.first + gt[5] * point.second + gt[4] * 0.5 + gt[5] * 0.5 + gt[3]
    );
}

std::vector<GeographicCoordinate> rowcol_to_projection(const GeoTransform& gt, const std::vector<ArrayCoordinate>& points) {
    std::vector<GeographicCoordinate> output{};
    output.reserve(points.size());
    
    std::transform(points.begin(), points.end(), std::back_inserter(output), [&gt](const auto& elem) {
        return rowcol_to_projection(gt, elem);
    });

    return output;
}

size_t rowcol_to_linear_index(const ArrayCoordinate& rowcol, const RasterShape& rs) {
    const size_t x_block = std::floor(rowcol.first / rs.blockxsize);
    const size_t y_block = std::floor(rowcol.second / rs.blockysize);
    const size_t linear_index = y_block * rs.size_x + x_block;

    return linear_index;
}

ArrayCoordinate linear_index_to_rowcol(const size_t linear_index, const RasterShape& rs) {
    const ArrayIndex_t row = std::floor(static_cast<double>(linear_index) / static_cast<double>(rs.size_x));
    const ArrayIndex_t col = linear_index - row*rs.size_x;

    return std::make_pair(col, row);
}

std::vector<GeographicCoordinate> read_points_from_geofile(OGRLayer* layer) {
    std::vector<GeographicCoordinate> points{};
    for (auto& feature : layer) {
        OGRGeometry* layer_geom = feature->GetGeometryRef();
        if (!layer_geom || !(wkbFlatten(layer_geom->getGeometryType()) == wkbPoint)) {
            std::cerr << "Detected feature with invalid/non-point geometry! Skipping...\n";
            continue;
        }

        OGRPoint* point = layer_geom->toPoint();
        points.push_back(std::make_pair((point->getX()), point->getY()));
    }

    return points;
}

CPLErr write_points_to_geofile(const std::string& path, const std::vector<GeographicCoordinate>& points, const std::vector<double>& values, const OGRSpatialReference* sref, const std::string output_driver, const std::string output_field_name = "sampled") {
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(output_driver.c_str());
    if (!driver) {
        std::cerr << "Could not use " << output_driver << " driver for writing\n";
        return CE_Failure;
    }

    GDALDataset* ds = driver->Create(path.c_str(), 0, 0, 0, GDT_Unknown, NULL);
    if (!ds) {
        std::cerr << "Could not create dataset using " << output_driver << " driver\n";
        return CE_Failure;
    }

    OGRSpatialReference* sref_output = sref->Clone();
    OGRLayer* layer = ds->CreateLayer("output", sref_output, wkbPoint, NULL);
    if (!layer) {
        std::cerr << "Could not create layer using " << output_driver << " driver\n";
        return CE_Failure;
    }

    OGRFieldDefn field(output_field_name.c_str(), OFTReal);
    field.SetWidth(10);
    if (layer->CreateField(&field) != OGRERR_NONE) {
        std::cerr << "Could not create field on output\n";
        return CE_Failure;
    }

    for (size_t i = 0; i < points.size(); i++) {
        OGRFeature* feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
        feature->SetField(output_field_name.c_str(), values[i]);
        
        OGRPoint pt(points[i].first, points[i].second);
        feature->SetGeometry(&pt);

        if (layer->CreateFeature(feature) != OGRERR_NONE) {
            std::cerr << "Failed to create feature!\n";
        }

        OGRFeature::DestroyFeature(feature);
    }

    GDALClose(ds);

    return CE_None;
}

void print_usage() {
    std::cout << "Usage: parsam [-l,--layer layername] [-b,--band band_number]\n";
    std::cout << "              [-cs,--cachesize size] [-of,--output_format format]\n";
    std::cout << "              input_raster input_points output_points\n\n";
    std::cout << "       -l, --layer: name of input point layer to sample. Defaults to first layer\n";
    std::cout << "       -b, --band: band of input raster to sample. Defaults to first band\n";
    std::cout << "       -cs, --cachesize: Size in megabyte for GDAL to use internally for caching\n";
    std::cout << "       -of, --output_format: Output format to use for result file. Needs to be a valid choice for GDAL tools\n\n";
    std::cout << "       input_raster: path pointing to input raster to sample points from\n";
    std::cout << "       input_points: path pointing to input vector layer containing points to be sampled\n";
    std::cout << "       output_points: path pointing to output vector layer containing sampled points\n";
}

void parse_cli_args(const int argc, const char* argv[], Arguments& args) {
    bool invalid_arg_found = false;
    size_t current_arg_idx = 1;
    size_t io_arg_idx = 0;
    std::array<std::string, 3> io_args{};
    while (current_arg_idx < argc) {
        const std::string current_arg{argv[current_arg_idx]};
        if (current_arg.compare("-l") == 0 || current_arg.compare("--layer") == 0) {
            args.point_layer = std::string{argv[current_arg_idx + 1]};
            current_arg_idx += 2;
            continue;
        }
        
        if (current_arg.compare("-b") == 0 || current_arg.compare("--band") == 0) {
            args.raster_band = std::stoi(argv[current_arg_idx + 1]);
            current_arg_idx += 2;
            continue;
        }

        if (current_arg.compare("-cs") == 0 || current_arg.compare("--cachesize") == 0) {
            args.cache_size = std::stoi(argv[current_arg_idx + 1]);
            current_arg_idx += 2;
            continue;
        }

        if (current_arg.compare("-of") == 0 || current_arg.compare("--output_format") == 0) {
            args.output_format = std::string{argv[current_arg_idx + 1]};
            current_arg_idx += 2;
            continue;
        }

        if (current_arg.compare(0, 1, "-") == 0 || current_arg.compare(0, 2, "--") == 0) {
            std::cerr << "Unknown argument: " << current_arg << "\n";
            invalid_arg_found = true;
            break;
        }
        
        io_args[io_arg_idx++] = current_arg;
        current_arg_idx++;
    }

    if (!invalid_arg_found) {
        args.input_raster = io_args[0];
        args.input_points = io_args[1];
        args.output_points = io_args[2];
    }
}

bool validate_parsed_arguments(const Arguments& args) {
    const bool raster_band_is_positive = args.raster_band > 0;
    const bool cachesize_is_non_negative = args.cache_size >= 0;
    const bool input_points_are_present = !args.input_points.empty();
    const bool input_raster_is_present = !args.input_raster.empty();
    const bool output_points_are_present = !args.output_points.empty();

    return raster_band_is_positive &&
           cachesize_is_non_negative &&
           input_points_are_present &&
           input_raster_is_present &&
           output_points_are_present;
}

bool check_crs_match(const OGRSpatialReference* raster_sref, const OGRSpatialReference* point_sref) {
    if (!raster_sref || !point_sref) {
        std::cerr << "Could not read spatial reference from inputs\n";
        return false;
    }
    if (point_sref->IsSame(raster_sref) == 0) {
        std::cerr << "Spatial reference systems of inputs do not match\n";
        return false;
    }

    return true;
}

void sort_proxy_by_blocks(std::vector<size_t>& point_proxy, const std::vector<ArrayCoordinate>& points_rowcol, const RasterShape& rs) {
    std::sort(point_proxy.begin(), point_proxy.end(), [&](const auto& idx_a, const auto& idx_b) {
        const auto a = points_rowcol[idx_a];
        const size_t a_linear_index = rowcol_to_linear_index(a, rs);

        const auto b = points_rowcol[idx_b];
        const size_t b_linear_index = rowcol_to_linear_index(b, rs);

        // TODO: Add check for out-of-bounds points and always return larger than other point.
        //       This way, all out-of-bounds points should go to the end.
        return a_linear_index < b_linear_index;
    });
}

std::vector<std::pair<size_t, size_t>> calculate_nonempty_window_indices(const std::vector<size_t>& point_proxy, const std::vector<ArrayCoordinate>& points_rowcol, const RasterShape& rs) {
    size_t starting_boundary_index = 0;
    size_t starting_linear_index = rowcol_to_linear_index(points_rowcol[point_proxy[0]], rs);
    std::vector<std::pair<size_t, size_t>> window_indices{};
    
    // This marks the initial window to sample from
    window_indices.push_back(std::make_pair(starting_boundary_index, starting_linear_index));
    
    for (size_t i = 0; i < point_proxy.size(); i++) {
        const size_t current_linear_index = rowcol_to_linear_index(points_rowcol[point_proxy[i]], rs);

        if (starting_linear_index != current_linear_index) {
            starting_boundary_index = i;
            starting_linear_index = current_linear_index;

            window_indices.push_back(std::make_pair(starting_boundary_index, starting_linear_index));
        }
    }

    return window_indices;
}

RasterShape raster_shape_from_band(GDALRasterBand* band) {
    int blockxsize, blockysize;
    band->GetBlockSize(&blockxsize, &blockysize);
    
    return RasterShape {
        static_cast<size_t>(blockxsize),
        static_cast<size_t>(blockysize),
        static_cast<size_t>(band->GetXSize()),
        static_cast<size_t>(band->GetYSize()),
        static_cast<size_t>((band->GetXSize() + blockxsize - 1) / blockxsize),
        static_cast<size_t>((band->GetYSize() + blockysize - 1) / blockysize)
    };
}

std::vector<double> sample_points_from_band(GDALRasterBand* band, const std::vector<ArrayCoordinate>& points_rowcol, const std::vector<size_t>& points_proxy, const std::vector<std::pair<size_t, size_t>> window_indices, const RasterShape& rs) {
    const size_t max_point_index = points_rowcol.size();
    
    std::vector<double> points_sampled{};
    points_sampled.resize(points_rowcol.size());
    
    thread_local std::vector<float> block_buffer{};
    block_buffer.reserve(rs.blockxsize*rs.blockysize);
    
    for (const auto [proxy_index, linear_index] : window_indices) {
        size_t current_point = proxy_index;
        const auto [block_x, block_y] = linear_index_to_rowcol(linear_index, rs);
        int valid_blocksize_x{}, valid_blocksize_y{};
        band->GetActualBlockSize(block_x, block_y, &valid_blocksize_x, &valid_blocksize_y);

        if (band->ReadBlock(block_x, block_y, block_buffer.data()) != CE_None) {
            throw std::runtime_error("Error reading block!");
        }

        // The left edge is always multiples of the full blocksize from the origin
        // while the right edge may not be a full block size away due to partial blocks
        const size_t block_xmin = block_x*rs.blockxsize, block_xmax = block_xmin + valid_blocksize_x;
        const size_t block_ymin = block_y*rs.blockysize, block_ymax = block_ymin + valid_blocksize_y;
        
        while (true) {
            const auto point_index = points_proxy[current_point];
            const auto [glob_x, glob_y] = points_rowcol[point_index];
            const bool within_x_bounds = block_xmin <= glob_x && glob_x < block_xmax;
            const bool within_y_bounds = block_ymin <= glob_y && glob_y < block_ymax;

            if (!within_x_bounds || !within_y_bounds || current_point == max_point_index - 1) {
                break;
            }

            const size_t local_x = glob_x - block_xmin, local_y = glob_y - block_ymin;
            points_sampled[point_index] = block_buffer[local_y*rs.blockxsize + local_x];
            current_point++;
        }
    }

    return points_sampled;
}

int main(int argc, const char* argv[]) {
    if (argc <= 3) {
        print_usage();
        return 0;
    };

    Arguments args{};
    parse_cli_args(argc, argv, args);
    if (!validate_parsed_arguments(args)) {
        if (args.raster_band < 1) {
            std::cerr << "Raster band needs to be larger than or equal to 1\n";
        } 
        else if (args.cache_size < 0) {
            std::cerr << "Cache size needs to be larger than or equal to zero\n";
        }
        else {
            std::cerr << "Missing ";
            if (args.input_points.empty()) std::cerr << "input points, ";
            if (args.input_raster.empty()) std::cerr << "input raster, ";
            if (args.output_points.empty()) std::cerr << "output points";
            std::cerr << "\n";
        }

        print_usage();
        return 1;
    }

    GDALAllRegister();

    // Limit GDAL's internal caching size. At 64MB, there does not seem to be much benefit from more
    CPLSetConfigOption("GDAL_CACHEMAX", std::to_string(args.cache_size).c_str());

    const GDALAccess eAccess = GA_ReadOnly;
    GDALDatasetUniquePtr input_raster = GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpen(args.input_raster.c_str(), eAccess)));
    if(!input_raster) {
        std::cerr << "Could not open input raster at " << args.input_raster << "\n";
        return 1;
    }

    GDALDatasetUniquePtr input_points = GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpenEx(args.input_points.c_str(), GA_ReadOnly | GDAL_OF_VECTOR, NULL, NULL, NULL)));
    if (!input_points) {
        std::cerr << "Could not open input points at " << args.input_points << "\n";
        return 1;
    }

    GeoTransform gt_ir{};
    if (input_raster->GetGeoTransform(gt_ir.data()) != CE_None) {
        std::cerr << "Could not read geotransform from input raster\n";
        return 1;
    }

    GDALRasterBand* band = input_raster->GetRasterBand(args.raster_band);
    if (!band) {
        std::cerr << "Could not read band " << args.raster_band << " from raster\n";
        return 1;
    }

    OGRLayer* point_layer{};
    if (!args.point_layer.empty()) {
        point_layer = input_points->GetLayerByName(args.point_layer.c_str());
    } else {
        point_layer = input_points->GetLayer(0);
    }

    if (!point_layer) {
        std::cerr << "Could not read layer from input point dataset\n";
        return 1;
    }

    if (!check_crs_match(input_raster->GetSpatialRef(), point_layer->GetSpatialRef())) {
        return 1;
    }

    const RasterShape rs = raster_shape_from_band(band);

    const auto points = read_points_from_geofile(point_layer);
    const auto points_rowcol = projection_to_rowcol(gt_ir, points);

    std::vector<size_t> points_proxy{};
    points_proxy.resize(points_rowcol.size());
    std::iota(points_proxy.begin(), points_proxy.end(), 0);
    sort_proxy_by_blocks(points_proxy, points_rowcol, rs);

    const auto window_indices = calculate_nonempty_window_indices(points_proxy, points_rowcol, rs);

    const auto points_sampled = sample_points_from_band(band, points_rowcol, points_proxy, window_indices, rs);

    CPLErr err = write_points_to_geofile(args.output_points, points, points_sampled, point_layer->GetSpatialRef(), args.output_format);
    if (err != CE_None) {
        std::cerr << "Error writing result\n";
        return 1;
    }

    return 0;
}