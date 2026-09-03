#!/usr/bin/env bash
set -euo pipefail

profile=${1:-test}
mode=${2:-run}
engine=${TEST_ENGINE:-docker}
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)

profiles=(test test-sb test-tm test-rustar test-bqc-star test-annevo)

# "<nextflow -profile>|<extra CLI args>". test-annevo reuses the `test` profile and
# only overrides what the ANNEVO path needs, so the two can never drift apart.
resolve() {
  case "$1" in
    test-annevo)
      echo "test|--annevo_annotation $root/assets/test/test_data/e2e/annevo.gff3 --output_dir $root/test_results/annevo"
      ;;
    *) echo "$1|" ;;
  esac
}

case "$profile" in
  all) selected=("${profiles[@]}") ;;
  *)
    if [[ " ${profiles[*]} " == *" $profile "* ]]; then
      selected=("$profile")
    else
      echo "usage: $0 {${profiles[*]// /|}|all} [run|verify]" >&2; exit 2
    fi
    ;;
esac

if [[ "$mode" != verify ]]; then
  nextflow_bin=${NEXTFLOW_BIN:-$(command -v nextflow || true)}
  if [[ -z "$nextflow_bin" ]]; then
    echo "Nextflow is not on PATH; run 'mamba activate nextflow' or set NEXTFLOW_BIN" >&2
    exit 2
  fi
  for name in "${selected[@]}"; do
    IFS='|' read -r base extra < <(resolve "$name")
    # shellcheck disable=SC2086
    "$nextflow_bin" run "$root" -profile "$base,$engine" -ansi-log false $extra
  done
fi

python "$root/assets/ci/e2e/verify.py" "$profile"
