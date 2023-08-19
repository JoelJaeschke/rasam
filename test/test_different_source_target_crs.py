import os
import tempfile
import subprocess
import filecmp

ASSETS_PATH = os.path.join(os.path.dirname(__file__), "assets")

def test_can_transform_different_input_crs_points(globals):
    raster_path = os.path.join(ASSETS_PATH, "input_raster_4326.tif")
    point_path = os.path.join(ASSETS_PATH, "sample_points.gpkg")
    reference_path = os.path.join(ASSETS_PATH, "reference_4326_from_3857_band1.csv")

    with tempfile.TemporaryDirectory() as tmp_path:
        output_path = os.path.join(tmp_path, "output.csv")
        result = subprocess.run([globals.get("UTIL_PATH"), raster_path, point_path, output_path, "-of", "CSV", "-l", "output_3857", "-s_srs", "EPSG:4326"], capture_output=True)

        assert filecmp.cmp(reference_path, output_path), "Files are not equal!"

    assert result.returncode == 0, "Tool should properly sample the input points"

def test_can_project_points_to_output_crs(globals):
    raster_path = os.path.join(ASSETS_PATH, "input_raster_4326.tif")
    point_path = os.path.join(ASSETS_PATH, "sample_points.gpkg")
    reference_path = os.path.join(ASSETS_PATH, "reference_3857_from_4326.csv")

    with tempfile.TemporaryDirectory() as tmp_path:
        output_path = os.path.join(tmp_path, "output.csv")
        result = subprocess.run([globals.get("UTIL_PATH"), raster_path, point_path, output_path, "-of", "CSV", "-l", "output_3857", "-s_srs", "EPSG:4326", "-t_srs", "EPSG:3857"], capture_output=True)

        assert filecmp.cmp(reference_path, output_path), "Files are not equal!"

    assert result.returncode == 0, "Tool should properly sample the input points"