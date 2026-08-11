#!/usr/bin/env bash
# EJ-230 center validation run — t0minidaq (24 threads).
exec "$(dirname "$0")/run_center.sh" --threads 24 "$@"
