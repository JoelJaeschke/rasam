# RASAM
A tool to sample raster files using vector inputs. This tool originated from necessity and the best alternative was `gdallocationinfo` which is barely useful for sampling large numbers of points. Using `rasterio` integrated sampling functionality is an alternative to this tool and there is still a need to do benchmarks of these two tools in order to figure out whether this tool has any benefits at all!

## How to use?
Usage is simple and all options are explained in the tool help. To get started, run
```shell
rasam <input_raster.tif> <input_points.tif> <output_file.{gpkg,csv,shp,...}>
```
The tool will automatically infer the output file type (if possible) and return all coordinates along with the sampled raster value.

## How to install?
You can build the tool manually using `cmake` by running
```shell
mkdir build
cd build && cmake ..
make
```
This will place the tool into `/build/src/rasam`. Alternatively, you will be able to install the tool using `conda` in the future!
