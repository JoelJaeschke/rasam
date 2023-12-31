#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gdal.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "cpl_error.h"
#include "ogr_feature.h"
#include "ogr_geometry.h"
#include "ogr_core.h"
#include "ogr_srs_api.h"
#include "ogr_spatialref.h"
#include "cpl_port.h"
#include "cpl_string.h"

using GeoTransform = std::array<double, 6>;
using ArrayIndex_t = size_t;
using Row = ArrayIndex_t;
using Column = ArrayIndex_t;
using ArrayCoordinate = std::pair<Column, Row>;
using GeographicIndex_t = double;
using Latitude = GeographicIndex_t;
using Longitude = GeographicIndex_t;
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
    std::string source_srs;
    std::string target_srs;

    Arguments(): raster_band(1),
                 cache_size(64),
                 point_layer(),
                 input_raster(),
                 input_points(),
                 output_points(),
                 output_format("GPKG"),
                 source_srs(),
                 target_srs() {};
    Arguments(const size_t raster_band,
              const size_t cache_size,
              const std::string point_layer,
              const std::string input_raster,
              const std::string input_points,
              const std::string output_points,
              const std::string output_format,
              const std::string source_srs,
              const std::string target_srs): raster_band(raster_band),
                                             cache_size(cache_size),
                                             point_layer(point_layer),
                                             input_raster(input_raster),
                                             input_points(input_points),
                                             output_points(output_points),
                                             output_format(output_format),
                                             source_srs(source_srs),
                                             target_srs(target_srs) {};
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

struct RasterBounds {
private:
    double ulx;
    double uly;
    double lrx;
    double lry;
public:
    RasterBounds(const GeoTransform& gt, const size_t rastersize_x, const size_t rastersize_y) {
        ulx = gt[0]; uly = gt[3] + static_cast<double>(rastersize_y)*gt[5];
        lrx = gt[0] + static_cast<double>(rastersize_x)*gt[1]; lry = gt[3];
    };
    ~RasterBounds() = default;

    std::array<GeographicIndex_t, 4> to_array() const {
        return {ulx, uly, lrx, lry};
    };
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

size_t rowcol_to_linear_index(const ArrayCoordinate& rowcol, const RasterShape& rs) {
    const size_t x_block = static_cast<size_t>(std::floor(rowcol.first / rs.blockxsize));
    const size_t y_block = static_cast<size_t>(std::floor(rowcol.second / rs.blockysize));
    const size_t linear_index = y_block * rs.size_x + x_block;

    return linear_index;
}

ArrayCoordinate linear_index_to_rowcol(const size_t linear_index, const RasterShape& rs) {
    const ArrayIndex_t row = static_cast<ArrayIndex_t>(std::floor(static_cast<double>(linear_index) / static_cast<double>(rs.size_x)));
    const ArrayIndex_t col = linear_index - row*rs.size_x;

    return std::make_pair(col, row);
}

std::vector<GeographicCoordinate> read_points_from_geofile(OGRLayer* layer, OGRCoordinateTransformation* transform) {
    OGRPoint* point = nullptr;
    OGRGeometry* layer_geom = nullptr;
    std::vector<GeographicCoordinate> points{};

    for (auto& feature : layer) {
        layer_geom = feature->GetGeometryRef();
        if (!layer_geom || !(wkbFlatten(layer_geom->getGeometryType()) == wkbPoint)) {
            std::cerr << "Detected feature with invalid/non-point geometry! Skipping...\n";
            continue;
        }

        point = layer_geom->toPoint();
        if (transform && point->transform(transform) != OGRERR_NONE) {
            std::cerr << "Failed to project point from point SRS to raster SRS!\n";
            continue;
        }
        points.push_back(std::make_pair(point->getX(), point->getY()));
    }

    return points;
}

OGRFieldType ogr_type_from_gdal_type(const GDALDataType dt) {
    if (dt == GDT_Float32 || dt == GDT_Float64) { return OFTReal; }
    if (dt == GDT_Int8 || dt == GDT_Int16 || dt == GDT_Int32 || dt == GDT_Int64 || 
        dt == GDT_UInt16 || dt == GDT_UInt32 || dt == GDT_UInt64) { return OFTInteger64; }
    
    return OFTMaxType;
}

template<typename T>
CPLErr write_points_to_geofile(const std::string& path, 
                               const std::vector<GeographicCoordinate>& points,
                               const std::vector<T>& values,
                               const GDALDataType raster_dt,
                               const OGRSpatialReference* sref,
                               OGRCoordinateTransformation* transform,
                               const std::string output_driver,
                               const std::string output_field_name = "sampled") {
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

    CPLStringList options{};
    if (output_driver.compare("CSV") == 0) {
        options.SetNameValue("GEOMETRY", "AS_XY");
    }

    OGRSpatialReference* sref_output = sref->Clone();
    OGRLayer* layer = ds->CreateLayer("output", sref_output, wkbPoint, options);
    if (!layer) {
        std::cerr << "Could not create layer using " << output_driver << " driver\n";
        return CE_Failure;
    }

    const OGRFieldType field_type = ogr_type_from_gdal_type(raster_dt);
    if (field_type == OFTMaxType) {
        std::cerr << "Could not determine field type!\n";
        return CE_Failure;
    }

    OGRFieldDefn field(output_field_name.c_str(), field_type);
    field.SetWidth(32);
    field.SetPrecision(13);

    if (layer->CreateField(&field) != OGRERR_NONE) {
        std::cerr << "Could not create field on output\n";
        return CE_Failure;
    }

    for (size_t i = 0; i < points.size(); i++) {
        OGRFeature* feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
        if constexpr (std::is_floating_point_v<T>) {
            feature->SetField(output_field_name.c_str(), static_cast<double>(values[i]));
        } else {
            feature->SetField(output_field_name.c_str(), static_cast<int>(values[i]));
        }

        OGRPoint pt{points[i].first, points[i].second};
        if (transform && pt.transform(transform) != OGRERR_NONE) {
            std::cerr << "Failed to project point from source SRS to specified target SRS!\n";
        };
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
    std::cout << "              [-s_srs,--source_srs srs] [-t_srs,--target_srs srs]\n";
    std::cout << "              input_raster input_points output_points\n\n";
    std::cout << "       -l, --layer: Name of input point layer to sample. Defaults to first layer.\n";
    std::cout << "       -b, --band: Band of input raster to sample. Defaults to first band.\n";
    std::cout << "       -cs, --cachesize: Size in megabyte for GDAL to use internally for caching.\n";
    std::cout << "       -of, --output_format: Output format to use for result file. Needs to be a valid choice for GDAL tools.\n";
    std::cout << "       -s_srs,--source_srs: Source spatial reference system to use for inputs. Note that input points will be\n"
              << "                            dynamically reprojected if necessary, but an error will be thrown if the raster\n"
              << "                            is not using this reference system.\n";
    std::cout << "       -t_srs,--target_srs: Target spatial reference system to project output points into.\n";
    std::cout << "\n";
    std::cout << "       input_raster: Path pointing to input raster to sample points from.\n";
    std::cout << "       input_points: Path pointing to input vector layer containing points to be sampled.\n";
    std::cout << "       output_points: Path pointing to output vector layer containing sampled points.\n";
}

void parse_cli_args(const int argc, const char* argv[], Arguments& args) {
    bool invalid_arg_found = false;
    size_t current_arg_idx = 1;
    size_t io_arg_idx = 0;
    std::array<std::string, 3> io_args{};
    while (current_arg_idx < static_cast<size_t>(argc)) {
        const std::string current_arg{argv[current_arg_idx]};
        if (current_arg.compare("-l") == 0 || current_arg.compare("--layer") == 0) {
            args.point_layer = std::string{argv[current_arg_idx + 1]};
            current_arg_idx += 2;
            continue;
        }
        
        if (current_arg.compare("-b") == 0 || current_arg.compare("--band") == 0) {
            const int raster_band = std::stoi(argv[current_arg_idx + 1]);
            if (raster_band < 1) {
                args.raster_band = 0;
            } else {
                args.raster_band = static_cast<size_t>(raster_band);
            }
            current_arg_idx += 2;
            continue;
        }

        if (current_arg.compare("-cs") == 0 || current_arg.compare("--cachesize") == 0) {
            const int cache_size = std::stoi(argv[current_arg_idx + 1]);
            if (cache_size < 0) {
                std::cerr << "Negative cachesize is not allowed!\n";
                invalid_arg_found = true;
                break;
            }

            args.cache_size = static_cast<size_t>(cache_size);
            current_arg_idx += 2;
            continue;
        }

        if (current_arg.compare("-of") == 0 || current_arg.compare("--output_format") == 0) {
            args.output_format = std::string{argv[current_arg_idx + 1]};
            current_arg_idx += 2;
            continue;
        }

        if (current_arg.compare("-s_srs") == 0 || current_arg.compare("--source_srs") == 0) {
            args.source_srs = std::string{argv[current_arg_idx + 1]};
            current_arg_idx += 2;
            continue;
        }

        if (current_arg.compare("-t_srs") == 0 || current_arg.compare("--target_srs") == 0) {
            args.target_srs = std::string{argv[current_arg_idx + 1]};
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
        std::cerr << "Could not read spatial reference from inputs\n";
        return false;
    }
    if (!point_sref->IsSame(raster_sref)) {
        std::cerr << "Spatial reference systems of inputs do not match\n";
        return false;
    }

    return true;
}

OGRErr determine_source_srs(const Arguments& args, OGRSpatialReference*& source_sref, const OGRSpatialReference* raster_sref, const OGRSpatialReference* point_sref) {
    if (!args.source_srs.empty()) {
        OGRSpatialReference user_input_srs;
        if (user_input_srs.SetFromUserInput(args.source_srs.c_str()) != OGRERR_NONE) {
            std::cerr << "Did not understand source spatial reference definition!\n";
            return OGRERR_UNSUPPORTED_SRS;
        };

        // TODO: This only checks whether they use the same geographical CRS and both are projected/non-projected. However,
        //       the projections may still be different!
        if (!(raster_sref->IsSameGeogCS(&user_input_srs) && raster_sref->IsProjected() == user_input_srs.IsProjected())) {
            std::cerr << "SRS mismatch between input raster and source raster definition. Are they using the same CRS?\n";
            return OGRERR_UNSUPPORTED_SRS;
        }

        source_sref = user_input_srs.Clone();
    } else {
        if (!check_crs_match(raster_sref, point_sref)) {
            std::cerr << "Mismatch between input raster and input point SRS! Either specify a custom one using '-s_srs' to reproject points\n"
                        << "on the fly to the raster SRS or reproject both to the same SRS before.\n";
            return OGRERR_UNSUPPORTED_SRS;
        }

        source_sref = raster_sref->Clone();
    }

    source_sref->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    return OGRERR_NONE;
}

OGRErr determine_target_srs(const Arguments& args, OGRSpatialReference*& target_sref, const OGRSpatialReference* source_sref) {
    if (!args.target_srs.empty()) {
        OGRSpatialReference user_target_srs;
        if (user_target_srs.SetFromUserInput(args.target_srs.c_str()) != OGRERR_NONE) {
            std::cerr << "Did not understand target spatial reference definition!\n";
            return OGRERR_UNSUPPORTED_SRS;
        }

        target_sref = user_target_srs.Clone();
    } else {
        target_sref = source_sref->Clone();
    }

    target_sref->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    return OGRERR_NONE;
}

std::vector<bool> calculate_points_in_bounds_mask(const std::vector<GeographicCoordinate>& points, const RasterBounds& rb) {
    std::vector<bool> in_bounds{};
    in_bounds.reserve(points.size());

    const auto raster_bounds = rb.to_array();

    std::transform(points.begin(), points.end(), std::back_inserter(in_bounds), [&raster_bounds](const auto& point) {
        const bool in_x_bounds = point.first > raster_bounds[0] && point.first < raster_bounds[2];
        const bool in_y_bounds = point.second > raster_bounds[1] && point.second < raster_bounds[3];
        return in_x_bounds && in_y_bounds;
    });

    return in_bounds;
}

void sort_proxy_by_blocks(std::vector<size_t>& point_proxy, const std::vector<ArrayCoordinate>& points_rowcol, const std::vector<bool>& valid_mask, const RasterShape& rs) {
    std::sort(point_proxy.begin(), point_proxy.end(), [&](const auto& idx_a, const auto& idx_b) {
        const auto a = points_rowcol[idx_a];
        const size_t a_linear_index = valid_mask[idx_a] ? rowcol_to_linear_index(a, rs) : std::numeric_limits<size_t>::max();

        const auto b = points_rowcol[idx_b];
        const size_t b_linear_index = valid_mask[idx_b] ? rowcol_to_linear_index(b, rs) : std::numeric_limits<size_t>::max();

        return a_linear_index < b_linear_index;
    });
}

std::vector<std::pair<size_t, size_t>> calculate_nonempty_window_indices(
    const std::vector<size_t>& point_proxy,
    const std::vector<ArrayCoordinate>& points_rowcol,
    const std::vector<bool>& valid_mask,
    const RasterShape& rs
) {
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

            if (valid_mask[point_proxy[i]]) {
                window_indices.push_back(std::make_pair(starting_boundary_index, starting_linear_index));
            } else {
                // Push a window with a special value for the linear index. This will be used to abort once all valid points are covered.
                window_indices.push_back(std::make_pair(starting_boundary_index, std::numeric_limits<size_t>::max()));
                break;
            }
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

template<typename T>
std::vector<T> sample_points_from_band(GDALRasterBand* band,
                                       const std::vector<ArrayCoordinate>& points_rowcol,
                                       const std::vector<size_t>& points_proxy,
                                       const std::vector<std::pair<size_t, size_t>>& window_indices,
                                       const RasterShape& rs)
{
    const size_t max_point_index = points_rowcol.size();
    
    T sample_initializer_value{};
    if constexpr (std::is_floating_point<T>::value) {
        sample_initializer_value = std::numeric_limits<T>::quiet_NaN();
    } else {
        sample_initializer_value = std::numeric_limits<T>::max();
    }

    std::vector<T> points_sampled(points_rowcol.size(), sample_initializer_value);
    
    std::vector<T> block_buffer{};
    block_buffer.reserve(rs.blockxsize*rs.blockysize);

    for (const auto& [proxy_index, linear_index] : window_indices) {
        if (linear_index == std::numeric_limits<size_t>::max()) {
            break;
        }

        size_t current_point = proxy_index;
        const auto [block_x, block_y] = linear_index_to_rowcol(linear_index, rs);
        int valid_blocksize_x{}, valid_blocksize_y{};
        if (band->GetActualBlockSize(static_cast<int>(block_x), static_cast<int>(block_y), &valid_blocksize_x, &valid_blocksize_y) != CE_None) {
            throw std::runtime_error("Error reading true blocksize!");
        };

        if (band->ReadBlock(static_cast<int>(block_x), static_cast<int>(block_y), block_buffer.data()) != CE_None) {
            throw std::runtime_error("Error reading block!");
        }

        // The left edge is always multiples of the full blocksize from the origin
        // while the right edge may not be a full block size away due to partial blocks
        const size_t block_xmin = block_x*rs.blockxsize, block_xmax = block_xmin + static_cast<size_t>(valid_blocksize_x);
        const size_t block_ymin = block_y*rs.blockysize, block_ymax = block_ymin + static_cast<size_t>(valid_blocksize_y);
        
        while (true) {
            const auto point_index = points_proxy[current_point];
            const auto [glob_x, glob_y] = points_rowcol[point_index];
            const bool within_x_bounds = block_xmin <= glob_x && glob_x < block_xmax;
            const bool within_y_bounds = block_ymin <= glob_y && glob_y < block_ymax;

            if (!within_x_bounds || !within_y_bounds || current_point == max_point_index) {
                break;
            }

            const size_t local_x = glob_x - block_xmin, local_y = glob_y - block_ymin;
            points_sampled[point_index] = block_buffer[local_y*rs.blockxsize + local_x];
            current_point++;
        }
    }

    return points_sampled;
}

CPLErr dispatch_sample_and_write(const GDALDataType raster_dt,
                                 GDALRasterBand* band,
                                 const std::vector<GeographicCoordinate>& points,
                                 const std::vector<ArrayCoordinate>& points_rowcol,
                                 const std::vector<size_t>& points_proxy,
                                 const std::vector<std::pair<size_t, size_t>>& window_indices,
                                 const RasterShape& rs,
                                 const OGRSpatialReference* target_sref,
                                 OGRCoordinateTransformation* target_point_transform,
                                 const Arguments& args)
{
    switch (raster_dt) {
        case GDT_Float32:
            {
                const auto points_sampled = sample_points_from_band<float>(band, points_rowcol, points_proxy, window_indices, rs);
                return write_points_to_geofile<float>(args.output_points, points, points_sampled, raster_dt, target_sref, target_point_transform, args.output_format);
            }
        case GDT_Float64:
            {
                const auto points_sampled = sample_points_from_band<double>(band, points_rowcol, points_proxy, window_indices, rs);
                return write_points_to_geofile<double>(args.output_points, points, points_sampled, raster_dt, target_sref, target_point_transform, args.output_format);
            }
        case GDT_UInt16:
            {
                const auto points_sampled = sample_points_from_band<uint16_t>(band, points_rowcol, points_proxy, window_indices, rs);
                return write_points_to_geofile<uint16_t>(args.output_points, points, points_sampled, raster_dt, target_sref, target_point_transform, args.output_format);
            }
        case GDT_UInt32:
            {
                const auto points_sampled = sample_points_from_band<uint32_t>(band, points_rowcol, points_proxy, window_indices, rs);
                return write_points_to_geofile<uint32_t>(args.output_points, points, points_sampled, raster_dt, target_sref, target_point_transform, args.output_format);
            }
        case GDT_UInt64:
            {
                const auto points_sampled = sample_points_from_band<uint64_t>(band, points_rowcol, points_proxy, window_indices, rs);
                return write_points_to_geofile<uint64_t>(args.output_points, points, points_sampled, raster_dt, target_sref, target_point_transform, args.output_format);
            }
        case GDT_Int8:
            {
                const auto points_sampled = sample_points_from_band<int8_t>(band, points_rowcol, points_proxy, window_indices, rs);
                return write_points_to_geofile<int8_t>(args.output_points, points, points_sampled, raster_dt, target_sref, target_point_transform, args.output_format);
            }
        case GDT_Int16:
            {
                const auto points_sampled = sample_points_from_band<int16_t>(band, points_rowcol, points_proxy, window_indices, rs);
                return write_points_to_geofile<int16_t>(args.output_points, points, points_sampled, raster_dt, target_sref, target_point_transform, args.output_format);
            }
        case GDT_Int32:
            {
                const auto points_sampled = sample_points_from_band<int32_t>(band, points_rowcol, points_proxy, window_indices, rs);
                return write_points_to_geofile<int32_t>(args.output_points, points, points_sampled, raster_dt, target_sref, target_point_transform, args.output_format);
            }
        case GDT_Int64:
            {
                const auto points_sampled = sample_points_from_band<int64_t>(band, points_rowcol, points_proxy, window_indices, rs);
                return write_points_to_geofile<int64_t>(args.output_points, points, points_sampled, raster_dt, target_sref, target_point_transform, args.output_format);
            }
        default:
            return CE_Failure;
    }
}

int main(int argc, const char* argv[]) {
    if (argc <= 3 || std::string(argv[1]).compare("-h") == 0) {
        print_usage();
        return 0;
    };

    Arguments args{};
    parse_cli_args(argc, argv, args);
    if (!validate_parsed_arguments(args)) {
        if (args.raster_band < 1) {
            std::cerr << "Raster band needs to be larger than or equal to 1\n";
        } else {
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

    GDALDatasetUniquePtr input_raster = GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpen(args.input_raster.c_str(), GA_ReadOnly)));
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

    GDALRasterBand* band = input_raster->GetRasterBand(static_cast<int>(args.raster_band));
    if (!band) {
        std::cerr << "Could not read band " << args.raster_band << " from raster\n";
        return 1;
    }

    OGRLayer* point_layer = nullptr;
    if (!args.point_layer.empty()) {
        point_layer = input_points->GetLayerByName(args.point_layer.c_str());
    } else {
        point_layer = input_points->GetLayer(0);
    }

    if (!point_layer) {
        std::cerr << "Could not read layer from input point dataset\n";
        return 1;
    }

    // If the user supplied custom source and target spatial reference systems,
    // check them early so we can fail early if something mismatches
    const OGRSpatialReference* raster_sref = input_raster->GetSpatialRef();
    const OGRSpatialReference* point_sref = point_layer->GetSpatialRef();

    OGRSpatialReference* source_sref = nullptr;
    if (determine_source_srs(args, source_sref, raster_sref, point_sref) != OGRERR_NONE) {
        return 1;
    }

    OGRSpatialReference* target_sref = nullptr;
    if (determine_target_srs(args, target_sref, source_sref) != OGRERR_NONE) {
        return 1;
    }

    // Transform input points during reading if they are using a different CRS. Just pass
    // a matching transform to the read_points_from_geofile function
    OGRCoordinateTransformation* source_point_transform = nullptr;
    if (!point_sref->IsSame(source_sref)) {
        source_point_transform = OGRCreateCoordinateTransformation(point_sref, source_sref);

        if (!source_point_transform) {
            std::cerr << "Failed to find a transform for going from point SRS to specified source SRS!\n";
            return 1;
        }
    }

    // Same applies for output SRS. Figure out early if we can even do this!
    OGRCoordinateTransformation *target_point_transform = nullptr;
    if (!source_sref->IsSame(target_sref)) {
        target_point_transform = OGRCreateCoordinateTransformation(source_sref, target_sref);

        if (!target_point_transform) {
            std::cerr << "Failed to find a transform for going from raster SRS to specified target SRS!\n";
            return 1;
        }
    }

    const RasterShape rs = raster_shape_from_band(band);
    const RasterBounds rb(
        gt_ir,
        static_cast<size_t>(input_raster->GetRasterXSize()),
        static_cast<size_t>(input_raster->GetRasterYSize())
    );

    const auto points = read_points_from_geofile(point_layer, source_point_transform);
    const auto in_bounds = calculate_points_in_bounds_mask(points, rb);
    const auto points_rowcol = projection_to_rowcol(gt_ir, points);

    std::vector<size_t> points_proxy(points_rowcol.size());
    std::iota(points_proxy.begin(), points_proxy.end(), 0);
    sort_proxy_by_blocks(points_proxy, points_rowcol, in_bounds, rs);

    const auto window_indices = calculate_nonempty_window_indices(points_proxy, points_rowcol, in_bounds, rs);

    if (dispatch_sample_and_write(band->GetRasterDataType(), band, points, points_rowcol, points_proxy,
                                  window_indices, rs, target_sref, target_point_transform, args) != CE_None)
    {
        std::cerr << "Error writing result\n";
        return 1;
    }

    return 0;
}