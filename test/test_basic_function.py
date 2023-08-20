import os
import tempfile
import subprocess
import filecmp

ASSETS_PATH = os.path.join(os.path.dirname(__file__), "assets")

def test_tool_shows_help_without_args(globals):
    result = subprocess.run([globals.get("UTIL_PATH")], capture_output=True)
    assert result.returncode == 0, "Tool should return success when called without any arguments"
    assert len(result.stdout) > 0, "Tool should display help text when called without arguments"

def test_tool_throws_error_with_invalid_args(globals):
    result = subprocess.run([globals.get("UTIL_PATH"), "input1", "input2", "--myarg", "invalid"], capture_output=True)
    assert result.returncode == 1, "Tool should throw error when called with invalid arguments"
    assert len(result.stderr) > 0, "Tool should display error messages"
    assert len(result.stdout) > 0, "Tool should display help/usage"

def test_can_sample_basic_raster(globals):
    raster_path = os.path.join(ASSETS_PATH, "input_raster_4326.tif")
    point_path = os.path.join(ASSETS_PATH, "sample_points.gpkg")
    reference_path = os.path.join(ASSETS_PATH, "reference_4326_band1.csv")

    with tempfile.TemporaryDirectory() as tmp_path:
        output_path = os.path.join(tmp_path, "output.csv")
        result = subprocess.run([globals.get("UTIL_PATH"), raster_path, point_path, output_path, "-of", "CSV"], capture_output=True)

        assert filecmp.cmp(reference_path, output_path), "Files are not equal!"

    assert result.returncode == 0, "Tool should properly sample the input points"

def test_can_sample_raster_with_different_layer(globals):
    raster_path = os.path.join(ASSETS_PATH, "input_raster_3857.tif")
    point_path = os.path.join(ASSETS_PATH, "sample_points.gpkg")
    reference_path = os.path.join(ASSETS_PATH, "reference_3857.csv")

    with tempfile.TemporaryDirectory() as tmp_path:
        output_path = os.path.join(tmp_path, "output.csv")    
        result = subprocess.run([globals.get("UTIL_PATH"), raster_path, point_path, output_path, "-of", "CSV", "-l", "output_3857"], capture_output=True)

        assert filecmp.cmp(output_path, reference_path), "Files are not equal!"

    assert result.returncode == 0, "Tool should properly sample the input points"

def test_can_sample_raster_with_different_band(globals):
    raster_path = os.path.join(ASSETS_PATH, "input_raster_4326.tif")
    point_path = os.path.join(ASSETS_PATH, "sample_points.gpkg")
    reference_path = os.path.join(ASSETS_PATH, "reference_4326_band2.csv")

    with tempfile.TemporaryDirectory() as tmp_path:
        output_path = os.path.join(tmp_path, "output.csv")    
        result = subprocess.run([globals.get("UTIL_PATH"), raster_path, point_path, output_path, "-of", "CSV", "-b", "2"], capture_output=True)

        assert filecmp.cmp(output_path, reference_path), "Files are not equal!"
    
    assert result.returncode == 0, "Tool should properly sample the input points"

def test_can_sample_raster_with_different_tiling(globals):
    raster_path = os.path.join(ASSETS_PATH, "input_raster_3857_tiled.tif")
    point_path = os.path.join(ASSETS_PATH, "sample_points.gpkg")
    reference_path = os.path.join(ASSETS_PATH, "reference_3857.csv")

    with tempfile.TemporaryDirectory() as tmp_path:
        output_path = os.path.join(tmp_path, "output.csv")    
        result = subprocess.run([globals.get("UTIL_PATH"), raster_path, point_path, output_path, "-of", "CSV", "-l", "output_3857"], capture_output=True)

        assert filecmp.cmp(output_path, reference_path), "Files are not equal!"
    
    assert result.returncode == 0, "Tool should properly sample the input points"

def test_can_handle_integer_data(globals):
    raster_path = os.path.join(ASSETS_PATH, "input_raster_4326_int.tif")
    point_path = os.path.join(ASSETS_PATH, "sample_points.gpkg")
    reference_path = os.path.join(ASSETS_PATH, "reference_4326_band1_asint.csv")

    with tempfile.TemporaryDirectory() as tmp_path:
        output_path = os.path.join(tmp_path, "output.csv")
        result = subprocess.run([globals.get("UTIL_PATH"), raster_path, point_path, output_path, "-of", "CSV"], capture_output=True)

        assert filecmp.cmp(reference_path, output_path), "Files are not equal!"

    assert result.returncode == 0, "Tool should properly sample the input points"