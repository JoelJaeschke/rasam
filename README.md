# RASAM
A tool to sample raster files using vector inputs. This tool originated from necessity and the best alternative was `gdallocationinfo` which is barely useful for sampling large numbers of points. Using `rasterio` integrated sampling functionality is an alternative to this tool and there is still a need to do benchmarks of these two tools in order to figure out whether this tool has any benefits at all!

## How to install?
The simplest way to install `rasam` is by using conda. Simply activate the environment you want to install the tool into and run:
```shell
conda install -c conda-forge rasam
```
Now, you should be able to run `rasam -h` and get the help displayed.

Alternatively, you can build the tool manually by cloning the repository and using `cmake` by running
```shell
mkdir build
cd build && cmake ..
make
```
This will place the tool into `/build/src/rasam`.

## How to use?
Usage is simple and all options are explained in the tool help. To get started, run
```shell
rasam <input_raster.tif> <input_points.{gpkg, csv, shp...}> <output_file.{gpkg,csv,shp,...}>
```
The tool will automatically infer the output file type (if possible) and return all coordinates along with the sampled raster value.
