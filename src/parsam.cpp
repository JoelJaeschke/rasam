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

struct Arguments {
public:
    size_t raster_band;
    size_t cache_size;
    std::string point_layer;
    std::string input_raster;
    std::string input_points;
    std::string output_points;

    Arguments(): raster_band(1),
                 cache_size(64),
                 point_layer(""),
                 input_raster(""),
                 input_points(""),
                 output_points("") {};
    Arguments(const size_t raster_band,
              const size_t cache_size,
              const std::string point_layer,
              const std::string input_raster,
              const std::string input_points,
              const std::string output_points): raster_band(raster_band),
                                                cache_size(cache_size),
                                                point_layer(point_layer),
                                                input_raster(input_raster),
                                                input_points(input_points),
                                                output_points(output_points) {};
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

std::pair<size_t, size_t> projection_to_rowcol(const GeoTransform& gt, const std::pair<double, double>& point) {
    return std::make_pair(
        static_cast<size_t>((point.first - gt[0]) / gt[1]),
        static_cast<size_t>((point.second - gt[3]) / gt[5])
    );
}

std::vector<std::pair<size_t, size_t>> projection_to_rowcol(const GeoTransform& gt, const std::vector<std::pair<double, double>>& points) {
    std::vector<std::pair<size_t, size_t>> output{};
    output.reserve(points.size());
    
    std::transform(points.begin(), points.end(), std::back_inserter(output), [&gt](const auto& elem) {
        return projection_to_rowcol(gt, elem);
    });

    return output;
}

std::pair<double, double> rowcol_to_projection(const GeoTransform& gt, const std::pair<size_t, size_t>& point) {
    return std::make_pair(
        gt[1] * point.first + gt[2] * point.second + gt[1] * 0.5 + gt[2] * 0.5 + gt[0],
        gt[4] * point.first + gt[5] * point.second + gt[4] * 0.5 + gt[5] * 0.5 + gt[3]
    );
}

std::vector<std::pair<double, double>> rowcol_to_projection(const GeoTransform& gt, const std::vector<std::pair<size_t, size_t>>& points) {
    std::vector<std::pair<double, double>> output{};
    output.reserve(points.size());
    
    std::transform(points.begin(), points.end(), std::back_inserter(output), [&gt](const auto& elem) {
        return rowcol_to_projection(gt, elem);
    });

    return output;
}

std::vector<std::pair<double, double>> read_points_from_geofile(OGRLayer* layer) {
    std::vector<std::pair<double, double>> points{};
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

CPLErr write_points_to_geofile(const std::string& path, const std::vector<std::pair<double, double>>& points, const std::vector<double>& values, const OGRSpatialReference* sref, const std::string output_driver = "FlatGeobuf", const std::string output_field_name = "sampled") {
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
    std::cout << "Usage: parsam [-l,--layer layername] [-b,--band band_number] [-cs,--cachesize size] input_raster input_points output_points\n";
    std::cout << "\t-l,--layer: name of input point layer to sample. Defaults to first layer\n";
    std::cout << "\t-b,--band: band of input raster to sample. Defaults to first band\n";
    std::cout << "\t-cs,--cachesize: Size in megabyte for GDAL to use internally for caching\n";
    std::cout << "\tinput_raster: path pointing to input raster to sample points from\n";
    std::cout << "\tinput_points: path pointing to input vector layer containing points to be sampled\n";
    std::cout << "\toutput_points: path pointing to output vector layer containing sampled points\n";

    return;
}

void parse_cli_args(const int argc, const char* argv[], Arguments& args) {
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
        
        io_args[io_arg_idx++] = current_arg;
        current_arg_idx++;
    }

    args.input_raster = io_args[0];
    args.input_points = io_args[1];
    args.output_points = io_args[2];

    return;
}

bool validate_parsed_arguments(const Arguments& args) {
    const bool raster_band_is_positive = args.raster_band > 0;
    const bool input_points_are_present = !args.input_points.empty();
    const bool input_raster_is_present = !args.input_raster.empty();
    const bool output_points_are_present = !args.output_points.empty();

    return raster_band_is_positive && 
           input_points_are_present && 
           input_raster_is_present && 
           output_points_are_present;
}

bool check_crs_match(const OGRSpatialReference* raster_sref, const OGRSpatialReference* point_sref) {
    if (!raster_sref || !point_sref) {
        std::cerr << "Could not spatial reference from inputs\n";
        return false;
    }
    if (point_sref->IsSame(raster_sref) == 0) {
        std::cerr << "Spatial reference systems of inputs do not match\n";
        return false;
    }

    return true;
}

size_t rowcol_to_linear_index(const std::pair<size_t, size_t>& rowcol, const RasterShape& rs) {
    const size_t x_block = std::floor(rowcol.first / rs.blockxsize);
    const size_t y_block = std::floor(rowcol.second / rs.blockysize);
    const size_t linear_index = y_block * rs.size_x + x_block;

    return linear_index;
}

std::pair<size_t, size_t> linear_index_to_rowcol(const size_t linear_index, const RasterShape& rs) {
    const size_t row = std::floor(static_cast<double>(linear_index) / static_cast<double>(rs.size_x));
    const size_t col = linear_index - row*rs.size_x;

    return std::make_pair(row, col);
}

void sort_proxy_by_blocks(std::vector<size_t>& point_proxy, const std::vector<std::pair<size_t, size_t>>& points_rowcol, const RasterShape& rs) {
    std::sort(point_proxy.begin(), point_proxy.end(), [&](const auto& idx_a, const auto& idx_b) {
        const auto a = points_rowcol[idx_a];
        // const size_t a_x_block = std::floor(a.first / rs.blockxsize);
        // const size_t a_y_block = std::floor(a.second / rs.blockysize);
        // const size_t a_linear_index = a_y_block * rs.size_x + a_x_block;
        const size_t a_linear_index = rowcol_to_linear_index(a, rs);

        const auto b = points_rowcol[idx_b];
        // const size_t b_x_block = std::floor(b.first / rs.blockxsize);
        // const size_t b_y_block = std::floor(b.second / rs.blockysize);
        // const size_t b_linear_index = b_y_block * rs.size_x + b_x_block;
        const size_t b_linear_index = rowcol_to_linear_index(b, rs);

        // TODO: Add check for out-of-bounds points and always return larger than other point.
        //       This way, all out-of-bounds points should go to the end.
        return a_linear_index < b_linear_index;
    });
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

std::vector<double> sample_points_from_band(GDALRasterBand* band, const std::vector<std::pair<size_t, size_t>>& points_rowcol, const std::vector<size_t>& points_proxy, const RasterShape& rs) {
    size_t current_point{0};
    const size_t max_point_index = points_rowcol.size();
    
    std::vector<double> points_sampled{};
    points_sampled.resize(points_rowcol.size());
    
    std::vector<float> block_buffer{};
    block_buffer.reserve(rs.blockxsize*rs.blockysize);
    
    for (size_t block_y = 0; block_y < rs.num_blocks_y; block_y++) {
        for (size_t block_x = 0; block_x < rs.num_blocks_x; block_x++) {
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
            std::cerr << "Missing inputs detected. Please specify all of 'input_points', 'input_raster' and 'output_points'\n";
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

    RasterShape rs = raster_shape_from_band(band);

    const auto points = read_points_from_geofile(point_layer);
    const auto points_rowcol = projection_to_rowcol(gt_ir, points);

    std::vector<size_t> points_proxy{};
    points_proxy.resize(points_rowcol.size());
    std::iota(points_proxy.begin(), points_proxy.end(), 0);
    sort_proxy_by_blocks(points_proxy, points_rowcol, rs);

    const auto points_sampled = sample_points_from_band(band, points_rowcol, points_proxy, rs);

    CPLErr err = write_points_to_geofile(args.output_points, points, points_sampled, point_layer->GetSpatialRef());
    if (err != CE_None) {
        std::cerr << "Error writing result\n";
        return 1;
    }

    return 0;
}