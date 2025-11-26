#!/usr/bin/env bash

# author = "Alejandro Gonzales-Irribarren"
# email = "alejandrxgzi@gmail.com"
# github = "https://github.com/alejandrogzi"
# version: 0.0.7

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

BEAVER_ENV_FILE="$REPO_ROOT/modules/custom/beaver/run/environment.yml"
MAIN_ENV_FILE="$REPO_ROOT/assets/metassembly.yml"
PACKAGE_DIR="$REPO_ROOT/modules/custom/beaver/run/package"
BEAVER_CHANNEL_VALUE="file://$PACKAGE_DIR/"

CONDA_BASE_DIR="$REPO_ROOT/conda"
META_ENV_PREFIX="$CONDA_BASE_DIR/metassembly"

echo "---------------------------------------------------"
echo "> metassembly: Bulk RNA-seq pipeline for transcript assembly"
echo "> Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>"
echo "> Repository: https://github.com/alejandrogzi/metassembly"
echo "---------------------------------------------------"

echo
echo "> Configuring metassembly..."
echo

echo "Repo root:          $REPO_ROOT"
echo "Beaver env file:    $BEAVER_ENV_FILE"
echo "Main env file:      $MAIN_ENV_FILE"
echo "Beaver channel:     $BEAVER_CHANNEL_VALUE"
echo "Conda base dir:     $CONDA_BASE_DIR"
echo "Metassembly prefix: $META_ENV_PREFIX"
echo

# --- Helper: update first channel line in a YAML file --------
update_first_channel() {
  local file="$1"
  local new_value="$2"

  if [[ ! -f "$file" ]]; then
    echo "ERROR: cannot update channels in missing file:"
    echo "  $file"
    return 1
  fi

  local tmp_file
  tmp_file="$(mktemp)"

  awk -v newchan="$new_value" '
    BEGIN {
      in_channels = 0;
      done = 0;
    }

    /^channels:/ {
      in_channels = 1;
    }

    in_channels && /^\s*-\s/ && !done {
      # Capture leading whitespace
      match($0, /^(\s*)-\s/, m)
      indent = m[1]

      # Replace line with same indent but new channel
      printf "%s- %s\n", indent, newchan
      done = 1
      next
    }

    /^dependencies:/ {
      in_channels = 0
    }

    { print }
  ' "$file" > "$tmp_file"

  mv "$tmp_file" "$file"
}

# --- Helper: insert beaver channel at top of channels: in a YAML file ---
insert_beaver_channel_main() {
  local file="$1"
  local new_value="$2"

  if [[ ! -f "$file" ]]; then
    echo "ERROR: cannot update main env in missing file:"
    echo "  $file"
    return 1
  fi

  # If the channel is already present, do nothing
  if grep -Fq "$new_value" "$file"; then
    echo "Beaver channel already present in $file, skipping insert."
    return 0
  fi

  local tmp_file
  tmp_file="$(mktemp)"

  awk -v newchan="$new_value" '
    BEGIN {
      in_channels = 0;
      inserted = 0;
    }

    /^channels:/ {
      in_channels = 1;
    }

    in_channels && /^\s*-\s/ && !inserted {
      # Capture leading whitespace of the original first channel
      match($0, /^(\s*)-\s/, m)
      indent = m[1]

      # Insert the new channel first, with same indent
      printf "%s- %s\n", indent, newchan
      inserted = 1
      # Fall through and also print the original line
    }

    /^dependencies:/ {
      in_channels = 0
    }

    { print }
  ' "$file" > "$tmp_file"

  mv "$tmp_file" "$file"
}

# --- 1) Update beaver environment.yml first channel ----------------------
echo "Updating beaver environment.yml first channel..."
update_first_channel "$BEAVER_ENV_FILE" "$BEAVER_CHANNEL_VALUE"
echo "✔ Updated $BEAVER_ENV_FILE"
echo

# --- 2) Insert beaver channel into main metassembly.yml -----------------
# echo "Inserting beaver channel into main metassembly.yml..."
# insert_beaver_channel_main "$MAIN_ENV_FILE" "$BEAVER_CHANNEL_VALUE"
# echo "✔ Updated $MAIN_ENV_FILE"
# echo

# --- 3) Create work/conda directory -------------------------------------
echo "Creating conda base directory (if needed)..."
mkdir -p "$CONDA_BASE_DIR"
echo "✔ Ensured $CONDA_BASE_DIR exists"
echo

# --- 4) Determine which mamba to use ------------------------------------
MAMBA_CMD=""

if [[ -n "${MAMBA_EXE:-}" && -x "$MAMBA_EXE" ]]; then
  MAMBA_CMD="$MAMBA_EXE"
elif command -v mamba >/dev/null 2>&1; then
  MAMBA_CMD="$(command -v mamba)"
else
  echo "ERROR: Could not find mamba."
  echo " - MAMBA_EXE is not set to an executable,"
  echo " - and 'mamba' is not on PATH."
  echo
  echo "If you use fish, make sure you have something like:"
  echo "  set -gx MAMBA_EXE /path/to/mamba"
  echo "in your config, and that it's exported (-gx)."
  exit 1
fi

echo "Using mamba at: $MAMBA_CMD"
echo

# --- 5) Create metassembly conda env with mamba -------------------------
if [[ -d "$META_ENV_PREFIX" ]]; then
  echo "Conda env already exists at:"
  echo "  $META_ENV_PREFIX"
  echo "Skipping mamba create."
else
  echo "Creating metassembly conda env with mamba..."
  "$MAMBA_CMD" create -y \
    --prefix "$META_ENV_PREFIX" \
    -f "$MAIN_ENV_FILE"

  echo "✔ Created env at:"
  echo "  $META_ENV_PREFIX"
fi

echo
echo "Configure step finished."
