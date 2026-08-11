#!/usr/bin/env bash
# EJ-230 position scan — t0minidaq (24 threads).
exec "$(dirname "$0")/run_scan.sh" --threads 24 "$@"
