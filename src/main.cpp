#include <algorithm>
#include <cmath>
#include <cstddef>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <ogr_feature.h>
#include <ogr_geometry.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "gdal.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include <cpl_error.h>
#include <ogr_core.h>
#include <ogr_spatialref.h>

using GeoTransform = std::array<double, 6>;

// All the projection logic done according to: https://gis.stackexchange.com/questions/424356/gdal-geographic-coordinate-to-pixel-coordinate
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

std::vector<std::pair<double, double>> read_points_from_geofile(const std::string& path) {
    GDALDatasetUniquePtr ds = GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpenEx(path.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL)));
    if (!ds) {
        throw std::runtime_error("Could not open vector file containing points!");
    }

    OGRLayer* layer = ds->GetLayer(1);
    if (!layer) {
        throw std::runtime_error("Could not read layer from vector file!");
    }

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

void write_points_to_geofile(const std::string& path, const std::vector<std::pair<double, double>>& points, const std::vector<double>& values, const std::string output_driver = "GPKG", const std::string output_field_name = "sampled") {
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(output_driver.c_str());
    if (!driver) {
        throw std::runtime_error("Could not use " + output_driver + " driver for writing!");
    }

    GDALDataset* ds = driver->Create(path.c_str(), 0, 0, 0, GDT_Unknown, NULL);
    if (!ds) {
        throw std::runtime_error("Could not create dataset using " + output_driver + " driver");
    }

    OGRLayer* layer = ds->CreateLayer("output", NULL, wkbPoint, NULL);
    if (!layer) {
        throw std::runtime_error("Could not create layer using " + output_driver + " driver");
    }

    OGRFieldDefn field(output_field_name.c_str(), OFTReal);
    field.SetWidth(10);
    if (layer->CreateField(&field) != OGRERR_NONE) {
        throw std::runtime_error("Could not create field on output");
    }

    for (size_t i = 0; i < points.size(); i++) {
        OGRFeature* feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
        feature->SetField("sample", values[i]);
        
        OGRPoint pt(points[i].first, points[i].second);
        feature->SetGeometry(&pt);

        if (layer->CreateFeature(feature) != OGRERR_NONE) {
            std::cerr << "Failed to create feature!\n";
        }

        OGRFeature::DestroyFeature(feature);
    }

    GDALClose(ds);
}

int main(int argc, const char* argv[]) {
    if (argc != 4) {
        throw std::runtime_error("Not enough arguments! Usage: ./sample <raster_file> <points_file>");
    }

    const std::string raster_file = std::string(argv[1]);
    const std::string points_file = std::string(argv[2]);
    const std::string output_file = std::string(argv[3]);

    GDALAllRegister();
    const GDALAccess eAccess = GA_ReadOnly;
    GDALDatasetUniquePtr poDataset = GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpen(raster_file.c_str(), eAccess)));
    
    if(!poDataset) {
        throw std::runtime_error("Could not open raster file for reading!");
    }

    // Get geotransform from file
    GeoTransform gt{};
    if (poDataset->GetGeoTransform(gt.data()) == CE_None) {
        std::cout << "Origin of image: " << gt[0] << ", " << gt[3] << "\n";
        std::cout << "Pixel size: " << gt[1] << ", " << gt[5] << "\n";

        if (gt[5] < 0.0) {
            std::cout << "Image is north-up\n";
        } else {
            std::cout << "Image is south-up\n";
        }
    }

    // Read spatial reference
    OGRSpatialReference image_proj(poDataset->GetProjectionRef());
    OGRSpatialReference geo_proj;
    const OGRErr err = geo_proj.importFromEPSG(3857);
    if (err != OGRERR_NONE) {
        throw std::runtime_error("Unknown geographical projection given!");
    }

    // Project geographical coordinates to image coordinates and backward
    const std::pair<double, double> hamburg = std::make_pair(9.99142, 53.5456);
    const std::pair<size_t, size_t> hamburg_img = projection_to_rowcol(gt, hamburg);
    const std::pair<double, double> hamburg_geo = rowcol_to_projection(gt, hamburg_img);

    std::cout << "Hamburg in geographic coordinates: " << hamburg.first << ", " << hamburg.second << "\n";
    std::cout << "Hamburg in image coordinates: " << hamburg_img.first << ", " << hamburg_img.second << "\n";
    std::cout << "Hamburg reprojected to geo coordinates from image: " << hamburg_geo.first << ", " << hamburg_geo.second << "\n";
    
    // Read raster blocks and decide how many blocks exist
    GDALRasterBand* band = poDataset->GetRasterBand(1);
    
    const int size_x = band->GetXSize();
    const int size_y = band->GetXSize();
    int blocksize_x{}, blocksize_y{};
    band->GetBlockSize(&blocksize_x, &blocksize_y);

    int num_blocks_x = (band->GetXSize() + blocksize_x - 1) / blocksize_x;
    int num_blocks_y = (band->GetYSize() + blocksize_y - 1) / blocksize_y;

    std::cout << "Band has size: "<< band->GetXSize() << ", " << band->GetYSize() << "\n";
    std::cout << "Dataset contains blocks of shape: " << blocksize_x << ", " << blocksize_y << ".\n";
    std::cout << "There are " << num_blocks_x << "x" << num_blocks_y << " blocks in the dataset.\n";

    // Read number of valid pixels for partial blocks
    int valid_blocksize_x{}, valid_blocksize_y{};
    band->GetActualBlockSize(14, 0, &valid_blocksize_x, &valid_blocksize_y);
    std::cout << "Right-most block column has valid pixels: " << valid_blocksize_x << "x" << valid_blocksize_y << "\n";

    std::cout << "\n==============================================\n\n";

    // Read points
    //const auto points = read_points_from_csv(points_file);
    const auto points = read_points_from_geofile(points_file);
    const auto points_rowcol = projection_to_rowcol(gt, points);

    // Create one level of indirection. This is used as a proxy for the elements and is the thing that will be sorted
    // instead of sorting the point coordinates themselves. This way, we can structure the API such that we return
    // sampled points in the order of inputs
    std::vector<size_t> points_sort_proxy{};
    points_sort_proxy.resize(points_rowcol.size());
    std::iota(points_sort_proxy.begin(), points_sort_proxy.end(), 0);
    
    std::vector<double> points_sampled{};
    points_sampled.resize(points_rowcol.size());

    // Sort elements according to their linear index into their respective blocks and order by blocks.
    std::sort(points_sort_proxy.begin(), points_sort_proxy.end(), [&](const auto& idx_a, const auto& idx_b) {
        const auto a = points_rowcol[idx_a];
        const size_t a_x_block = std::floor(a.first / blocksize_x);
        const size_t a_y_block = std::floor(a.second / blocksize_y);
        const size_t a_linear_index = a_y_block * size_x + a_x_block;

        const auto b = points_rowcol[idx_b];
        const size_t b_x_block = std::floor(b.first / blocksize_x);
        const size_t b_y_block = std::floor(b.second / blocksize_y);
        const size_t b_linear_index = b_y_block * size_x + b_x_block;

        // TODO: Add check for out-of-bounds points and always return larger than other point.
        //       This way, all out-of-bounds points should go to the end.
        return a_linear_index < b_linear_index;
    });

    // Go through all blocks and iterate the points until one is out-of-bounds
    size_t current_point{0};
    const size_t max_point_index = points_rowcol.size();
    std::cout << "Max point index: " << max_point_index << "\n\n";

    std::vector<float> block_buffer{};
    block_buffer.reserve(blocksize_x*blocksize_y);
    for (size_t block_y = 0; block_y < num_blocks_y; block_y++) {
        for (size_t block_x = 0; block_x < num_blocks_x; block_x++) {
            int valid_blocksize_x{}, valid_blocksize_y{};
            band->GetActualBlockSize(block_x, block_y, &valid_blocksize_x, &valid_blocksize_y);

            CPLErr error = band->ReadBlock(block_x, block_y, block_buffer.data());
            if (error != CE_None) {
                throw std::runtime_error("Error reading block!");
            }

            // The left edge is always multiples of the full blocksize from the origin
            // while the right edge may not be a full block size away due to partial blocks
            const size_t block_xmin = block_x*blocksize_x, block_xmax = block_xmin + valid_blocksize_x;
            const size_t block_ymin = block_y*blocksize_y, block_ymax = block_ymin + valid_blocksize_y;
            
            while (true) {
                const auto point_index = points_sort_proxy[current_point];
                const auto [glob_x, glob_y] = points_rowcol[point_index];
                const bool within_x_bounds = block_xmin <= glob_x && glob_x < block_xmax;
                const bool within_y_bounds = block_ymin <= glob_y && glob_y < block_ymax;

                if (!within_x_bounds || !within_y_bounds || current_point == max_point_index - 1) {
                    break;
                }

                const size_t local_x = glob_x - block_xmin, local_y = glob_y - block_ymin;
                points_sampled[point_index] = block_buffer[local_y*blocksize_x + local_x];
                current_point++;

                // std::cout << "Global Bounding box: " << block_xmin << ", " << block_ymin << ", " << block_xmax << ", " << block_ymax << "\n";
                // std::cout << "Global point coordinate: " << glob_x << ", " << glob_y << "\n";
                // std::cout << "Local point coordinate: " << local_x << ", " << local_y << "\n";
            }
        }
    }

    // for (size_t i = 0; i < points_sampled.size(); i++) {
    //     const auto point_val = points_sampled[i];
    //     const auto point_rowcol = points_rowcol[i];
    //     const auto point_geo = rowcol_to_projection(gt, point_rowcol);
    //     std::cout << point_geo.first << "," << point_geo.second << "," << point_val << "\n";
    // }

    write_points_to_geofile(output_file, points, points_sampled);

    return 0;
}