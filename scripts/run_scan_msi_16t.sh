#!/usr/bin/env bash
# EJ-230 position scan — MSI (16 threads).
exec "$(dirname "$0")/run_scan.sh" --threads 16 "$@"
