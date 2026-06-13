#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NTHREADS=16
exec "$script_dir/run_scan.sh" --threads "$NTHREADS" --host msi "$@"
