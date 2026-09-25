# EXEC_31 — pre-integration commit classification

All commits have one unambiguous primary-purpose category; the two user-identified physics fixes remain physics-fix despite ancillary counter lines. No ambiguous commit gate is triggered. Main has not been edited.

| Commit | Category | Justification | Runtime cost |
|---|---|---|---|
| 8b64721 | experiment | Adds YPlus while TOP is present for D2; deliberate design comparison, not baseline correction. | Geometry change; do not integrate. |
| 2ddf205 | instrumentation | Adds the terminal boundary-state subset census without changing tracks. | Additional status lookup and map increments inside the terminal mutex. |
| f4d90a6 | physics-fix | Changes the air-wrap factory call to dielectric_metal with REFLECTIVITY=0.98. | Not diagnostic; physical trajectory lengths can change. |
| a7a9090 | instrumentation | Adds terminal identity/fate bookkeeping and observation of explicit kills, without altering those kills. | Per-track sets/maps, shared terminal mutex, kill-flag bookkeeping and export. |
| 218241a | physics-fix | Changes the world guard to preserve reflected photons; the added diagnostic counter is ancillary to this explicit physics fix. | Ancillary atomic counter; boundary-status lookup is required physics and must remain ON. |
| 0336ba9 | instrumentation | Records boundary outcomes by physical-volume pair without changing track decisions. | Per-boundary map lookup/increment under a shared mutex and volume-key construction. |
| 3f0808d | instrumentation | Adds the mutex-protected boundary census and RNG snapshots; no track-state or physics edits. | Map/mutex storage; reset/export and RNG-file I/O; per-encounter cost begins when wired by 0336ba9. |
| 391fa4e | experiment | Restores the original dielectric factory for D3; deliberate control, not production fix. | Physical control; preserve by tag, do not integrate. |
| 2e3872e | instrumentation | Extends current terminal-state observations to every fate, with no physical state writes. | Additional terminal maps under the existing shared mutex. |
| 413a8f0 | instrumentation | Centralizes shared diagnostic resets and summaries in master; leaves tracking physics unchanged. | Small lifecycle/I/O cost; removes repeated worker summaries. |

The EXEC_30 branch is listed separately because the worktree currently points there, while the requested log branch ends at D2. All diagnostic instrumentation can be excluded by an OFF-by-default CMake option, including its source files/actions and counter updates. Required physics boundary-status lookup stays unconditional. Ordinary event/scintillation totals and hit output remain part of the baseline.

The measured EXEC_30 V1 speedup was 1.468017 (3341.036238 s / 2275.883586 s, N=2000, seeds 26092601 8349041); source: build_exec30_20260910/audit/G0.json and cells/V1/cell.meta.json in the diagnostic worktree. Mutex contention is a plausible source-level explanation, not a separately measured profiler attribution. The prior report explicitly qualified this as an elapsed-time comparison, not an isolated scaling benchmark.
