#!/usr/bin/env bash
# EJ-230 center validation run — MSI (16 threads).
exec "$(dirname "$0")/run_center.sh" --threads 16 "$@"
