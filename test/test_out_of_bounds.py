import os
import tempfile
import subprocess
import filecmp

ASSETS_PATH = os.path.join(os.path.dirname(__file__), "assets")

def test_can_handle_out_of_bounds(globals):
    raster_path = os.path.join(ASSETS_PATH, "input_raster_4326.tif")
    point_path = os.path.join(ASSETS_PATH, "sample_points.gpkg")
    reference_path = os.path.join(ASSETS_PATH, "reference_4326_oob_band1.csv")

    with tempfile.TemporaryDirectory() as tmp_path:
        output_path = os.path.join(tmp_path, "output.csv")
        result = subprocess.run([globals.get("UTIL_PATH"), raster_path, point_path, output_path, "-of", "CSV", "-l", "sample_points_oob"], capture_output=True)

        assert filecmp.cmp(reference_path, output_path), "Files are not equal!"

    assert result.returncode == 0, "Tool should properly sample the input points"