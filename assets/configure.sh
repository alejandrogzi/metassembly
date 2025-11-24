#!/usr/bin/env bash

# author = "Alejandro Gonzales-Irribarren"
# email = "alejandrxgzi@gmail.com"
# github = "https://github.com/alejandrogzi"
# version: 0.0.7

set -euo pipefail

# Determine repo root as the parent directory of assets/
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ENV_FILE="$REPO_ROOT/modules/custom/beaver/run/environment.yml"
PACKAGE_DIR="$REPO_ROOT/modules/custom/beaver/run/package"
NEW_CHANNEL_VALUE="file://$PACKAGE_DIR/"

echo "Repo root:        $REPO_ROOT"
echo "Environment file: $ENV_FILE"
echo "Package dir:      $PACKAGE_DIR"
echo

# Sanity checks
if [[ ! -f "$ENV_FILE" ]]; then
  echo "ERROR: environment.yml not found at:"
  echo "  $ENV_FILE"
  exit 1
fi

# Process YAML and preserve indentation of the first channel
tmp_file="$(mktemp)"

awk -v newchan="$NEW_CHANNEL_VALUE" '
  BEGIN {
    in_channels = 0;
    done = 0;
  }

  /^channels:/ {
    in_channels = 1;
  }

  in_channels && /^\s*-\s/ && !done {
    # Capture leading whitespace
    match($0, /^(\s*)-/ , m)
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
' "$ENV_FILE" > "$tmp_file"

mv "$tmp_file" "$ENV_FILE"

echo "✔ Updated first channel to:"
echo "  - $NEW_CHANNEL_VALUE"
