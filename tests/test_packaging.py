# - Build a wheel in an isolated test destination.
import subprocess
from pathlib import Path
from zipfile import ZipFile

def test_distributed_notice(tmp_path):
    root = Path(__file__).parents[1]
    subprocess.run(["uv", "build", "--wheel", "--out-dir", str(tmp_path)], cwd=root,
                   check=True, capture_output=True)

# - Inspect distributable contents without importing source checkout.
    wheel = next(tmp_path.glob("*.whl"))
    with ZipFile(wheel) as archive:
        names = archive.namelist()
        notice = next(name for name in names if name.endswith("/licenses/UPSTREAM_LICENSE.txt"))
        metadata = next(name for name in names if name.endswith(".dist-info/METADATA"))

# - Verify upstream notice and license metadata survive packaging.
        assert archive.read(notice).decode() == (root / "UPSTREAM_LICENSE.txt").read_text()
        assert "Copyright (C) 2013 Stanford University" in archive.read(notice).decode()
        assert "License-File: UPSTREAM_LICENSE.txt" in archive.read(metadata).decode()
