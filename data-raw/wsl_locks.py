"""Export portable lock contents from a tested WSL environment as JSON.

Run with that environment's Python, passing the micromamba executable path.
Review the output before replacing the corresponding files in inst/wsl.
The three analysis tools retain their source-commit pins in the R installer.
"""

import importlib.metadata as metadata
import json
import re
import subprocess
import sys
import urllib.request

from packaging.tags import sys_tags
from packaging.utils import canonicalize_name, parse_wheel_filename


def export_locks():
    conda = subprocess.check_output([
        sys.argv[1], "list", "--prefix", sys.prefix,
        "--no-pip", "--explicit", "--sha256",
    ], text=True)
    # Keep the Python environments limited to Conda-forge packages.
    urls = [line for line in conda.splitlines()
            if line.startswith("https://conda.anaconda.org/conda-forge/")]
    assert urls and all(re.search(r"#[0-9a-f]{64}$", url) for url in urls)

    versions = {
        canonicalize_name(dist.metadata["Name"]): dist.version
        for dist in metadata.distributions()
        if canonicalize_name(dist.metadata["Name"])
        not in {"fastmixture", "adamixture", "clumppling"}
    }
    # Source-build backends are not retained by pip's isolated builds.
    if sys.version_info.minor == 11:
        versions["setuptools-scm"] = "8.3.1"
    elif sys.version_info.minor == 12:
        versions["poetry-core"] = "2.5.0"
    else:
        raise ValueError("Use the tested Python 3.11 or 3.12 environment")

    # Include compatible wheel hashes for Ubuntu 22.04 and newer (glibc 2.35+).
    tags = set()
    for tag in sys_tags():
        glibc = re.match(r"manylinux_2_(\d+)_x86_64$", tag.platform)
        if tag.platform == "any" or (glibc and int(glibc[1]) <= 35):
            tags.add(tag)
    requirements = []
    for name, version in sorted(versions.items()):
        with urllib.request.urlopen(
            f"https://pypi.org/pypi/{name}/{version}/json", timeout=60
        ) as response:
            release = json.load(response)
        wheels = []
        sources = []
        for item in release["urls"]:
            if item["yanked"]:
                continue
            if item["packagetype"] == "bdist_wheel":
                wheel_tags = parse_wheel_filename(item["filename"])[3]
                if tags.intersection(wheel_tags):
                    wheels.append(item["digests"]["sha256"])
            elif item["packagetype"] == "sdist":
                sources.append(item["digests"]["sha256"])
        hashes = sorted(set(wheels or sources))
        assert hashes, f"No compatible distribution for {name}=={version}"
        requirements.append(f"{name}=={version}" + "".join(
            f" --hash=sha256:{digest}" for digest in hashes
        ))
    print(json.dumps({
        "conda": "@EXPLICIT\n" + "\n".join(sorted(urls)) + "\n",
        "pip": "\n".join(requirements) + "\n",
    }))


if __name__ == "__main__":
    export_locks()
