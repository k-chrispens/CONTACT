# - Locate the public original CONTACT supplementary archive.
import argparse
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

# - Verify downloaded archive against pinned provenance.
URL = ("https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.2592/"
       "MediaObjects/41592_2013_BFnmeth2592_MOESM125_ESM.zip")
SHA256 = "9f3c5fe8ba0513d87d1d447d0622c64cf7119bb168cd19782f43400d11bbf594"

# - Extract only explicitly selected reference artifacts.
PREFIX = "new_DHFR_wt.updated_refine_001.pdb_modified.updated"
ARTIFACTS = {
    "dhfr.pdb": PREFIX + ".pdb", "pathways.txt": PREFIX + "_graph_0.30_100.txt",
    "maximal_pathways.txt": PREFIX + "_graph_0.30_100.txtRSDS_ONLY",
    "graph.log": "graph.log", "relief_by_residue.txt": "relief_by_residue.txt",
    "LICENSE.txt": "LICENSE.txt", "original_README": "README",
    "CreateGraph_stats.cpp": "CreateGraph_stats.cpp", "CreateGraph_stats.hpp": "CreateGraph_stats.hpp",
    "mainfile.cpp": "mainfile.cpp", "Makefile": "Makefile",
}

# - Refuse existing destinations and unsafe archive paths.
def fetch(destination):
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    with urlopen(URL, timeout=60) as response:
        payload = response.read(1_000_000)
    if hashlib.sha256(payload).hexdigest() != SHA256:
        raise ValueError("CONTACT supplementary archive checksum mismatch")
    destination.mkdir(parents=True, exist_ok=False)

# - Record source attribution and artifact checksums.
    checksums = {}
    with ZipFile(io.BytesIO(payload)) as archive:
        for filename, member in ARTIFACTS.items():
            data = archive.read("CONTACT/" + member)
            (destination / filename).write_bytes(data)
            checksums[filename] = hashlib.sha256(data).hexdigest()
    provenance = {"doi": "10.1038/nmeth.2592", "url": URL, "archive_sha256": SHA256,
                  "files": checksums, "notice": "Original MIT source; reference artifacts downloaded, not redistributed"}
    (destination / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return destination

# - Fetch only when explicitly invoked by users.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch original CONTACT DHFR references; never run bundled binaries")
    parser.add_argument("--output", required=True)
    print(fetch(parser.parse_args().output))
