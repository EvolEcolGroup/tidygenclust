#!/bin/sh
set -eu

# Arguments come from .tgc_wsl_install_tools(); locks live beside this script.
root=$1
fastmixture_hash=$2
clumppling_hash=$3
reset=$4
install_adamixture=$5
adamixture_hash=$6
lock_dir=$(dirname "$0")

# Check prerequisites before downloading anything or creating environments.
case "$(uname -r)" in
  *microsoft-standard*|*WSL2*|*wsl2*) ;;
  *) echo "Please use WSL2; convert your distribution with wsl --set-version <name> 2." >&2
     exit 1 ;;
esac
if [ "$(uname -m)" != "x86_64" ]; then
  echo "The Windows backend currently supports x86-64 computers, not Windows on ARM." >&2
  exit 1
fi
for tool in curl git sha256sum gcc g++ make tar gzip; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Missing $tool in WSL. In Ubuntu, run:" >&2
    echo "sudo apt update && sudo apt install -y build-essential curl git ca-certificates tar gzip" >&2
    exit 1
  fi
done

# Keep the manager, its cache and both analysis environments under one private
# Linux prefix, separate from any Windows or user-managed Conda installation.
bin_dir="$root/bin"
micromamba="$bin_dir/micromamba"
mamba_root="$root/micromamba"
core_env="$root/envs/ctidygenclust"
clumppling_env="$root/envs/cclumppling"

mkdir -p "$bin_dir" "$root/envs" "$mamba_root"
download_dir=$(mktemp -d)
trap 'rm -rf "$download_dir"' EXIT HUP INT TERM

# Check an existing executable too, so every installation uses this build.
mamba_version="2.9.0-0"
mamba_sha256="366cd9cd8be14df1ab8ed50352a82111082a36686b2d389fdb79a92c3fafb3e3"
if [ ! -x "$micromamba" ] ||
  ! printf '%s  %s\n' "$mamba_sha256" "$micromamba" |
    sha256sum --check --status 2>/dev/null; then
  echo "Installing micromamba in WSL..."
  # Pin the executable and verify it before running downloaded code.
  curl --fail --location --silent --show-error --retry 3 \
    "https://github.com/mamba-org/micromamba-releases/releases/download/$mamba_version/micromamba-linux-64" \
    --output "$download_dir/micromamba"
  printf '%s  %s\n' "$mamba_sha256" "$download_dir/micromamba" | sha256sum --check
  cp "$download_dir/micromamba" "$micromamba"
  chmod 755 "$micromamba"
fi

# Only an explicit reset removes environments; ordinary installs reuse them.
if [ "$reset" = "true" ]; then
  if [ -d "$core_env" ]; then
    "$micromamba" --root-prefix "$mamba_root" env remove \
      --yes --prefix "$core_env"
  fi
  if [ -d "$clumppling_env" ]; then
    "$micromamba" --root-prefix "$mamba_root" env remove \
      --yes --prefix "$clumppling_env"
  fi
fi

# The original tools need separate environments with different Python versions.
if [ ! -x "$core_env/bin/python" ]; then
  echo "Creating the ctidygenclust WSL environment..."
  "$micromamba" --root-prefix "$mamba_root" create \
    --yes --prefix "$core_env" \
    --file "$lock_dir/ctidygenclust-conda.lock"

else
  echo "The ctidygenclust WSL environment already exists; keeping it."
fi

if [ ! -x "$clumppling_env/bin/python" ]; then
  echo "Creating the cclumppling WSL environment..."
  "$micromamba" --root-prefix "$mamba_root" create \
    --yes --prefix "$clumppling_env" \
    --file "$lock_dir/cclumppling-conda.lock"

else
  echo "The cclumppling WSL environment already exists; keeping it."
fi

# Importing alone would silently keep an older revision. Check pip's source
# metadata too, so changing a requested commit actually updates the tool.
check_source='import importlib, importlib.metadata as m, json, sys
importlib.import_module(sys.argv[1])
source = json.loads(m.distribution(sys.argv[1]).read_text("direct_url.json") or "{}")
assert source.get("url") == sys.argv[2]
assert source.get("vcs_info", {}).get("requested_revision") == sys.argv[3]'

for tool in fastmixture adamixture clumppling; do
  pip_lock="$lock_dir/ctidygenclust-pip.lock"
  case "$tool" in
    fastmixture)
      env_path=$core_env
      repo="https://github.com/Rosemeis/fastmixture.git"
      revision=$fastmixture_hash ;;
    adamixture)
      [ "$install_adamixture" = "true" ] || continue
      env_path=$core_env
      repo="https://github.com/AI-sandbox/ADAMIXTURE.git"
      revision=$adamixture_hash ;;
    clumppling)
      env_path=$clumppling_env
      pip_lock="$lock_dir/cclumppling-pip.lock"
      repo="https://github.com/PopGenClustering/Clumppling.git"
      revision=$clumppling_hash ;;
  esac
  if ! "$micromamba" --root-prefix "$mamba_root" run --prefix "$env_path" \
    python -c "$check_source" "$tool" "$repo" "$revision" >/dev/null 2>&1; then
    echo "Installing $tool at revision $revision..."
    # Install the frozen runtime/build dependencies first. Do not let source
    # builds create a separate environment with newer, unpinned dependencies.
    "$micromamba" --root-prefix "$mamba_root" run --prefix "$env_path" \
      python -m pip install --no-cache-dir --require-hashes --no-deps \
      --no-build-isolation --requirement "$pip_lock"
    "$micromamba" --root-prefix "$mamba_root" run --prefix "$env_path" \
      python -m pip install --no-cache-dir --force-reinstall --no-deps \
      --no-build-isolation "git+$repo@$revision"
  fi
  "$micromamba" --root-prefix "$mamba_root" run --prefix "$env_path" \
    python -c "$check_source" "$tool" "$repo" "$revision"
done

# Also detect dependency conflicts in environments retained from earlier runs.
for env_path in "$core_env" "$clumppling_env"; do
  "$micromamba" --root-prefix "$mamba_root" run --prefix "$env_path" \
    python -m pip check
done

echo "tidygenclust WSL tools installed successfully."
