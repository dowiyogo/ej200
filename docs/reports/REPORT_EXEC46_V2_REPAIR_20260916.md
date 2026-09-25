# EXEC_46 v2 repair report

Date: 2026-09-16.

## 1. Diagnosis

The contaminated evidence is preserved at `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_BROKEN_20260916`. The 14 corrected cells had the same direct target:

`/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/sslg4`

The controlling code in `prepare_campaign.py` at commit `f801d76` was:

```python
runtime_sslg4 = (corrected_sslg4_source
                 if source_cell["material"] in ("EJ-200", "EJ-204")
                 and corrected_sslg4_source is not None else binary.parent / "sslg4")
(target / "sslg4").symlink_to(runtime_sslg4, target_is_directory=True)
```

`corrected_sslg4_source` is one function argument, and the loop never selected a source by position. The material condition only selected the same argument for both corrected materials; it did not distinguish EJ-200 from EJ-204 or any `cell_id`. Consequently, every corrected `target / "sslg4"` pointed to the same path. EJ-230 followed the `else` branch and pointed to `build_baseline/sslg4`.

The relevant commit diff was:

```diff
-        runtime_sslg4 = (corrected_sslg4_source
-                         if source_cell["material"] in ("EJ-200", "EJ-204")
-                         and corrected_sslg4_source is not None else binary.parent / "sslg4")
-        (target / "sslg4").symlink_to(runtime_sslg4, target_is_directory=True)
+        runtime_source = (CORRECTED_SSLG4_BY_MATERIAL[source_cell["material"]]
+                          if source_cell["material"] in CORRECTED_SSLG4_BY_MATERIAL
+                          and corrected_sslg4_source is not None
+                          else binary.parent / "sslg4")
+        runtime_link = runtime_links / source_cell["cell_id"]
+        runtime_link.symlink_to(runtime_source, target_is_directory=True)
+        runtime_sslg4 = runtime_link
+        local_sslg4 = target / "sslg4"
+        local_sslg4.symlink_to(runtime_sslg4, target_is_directory=True)
+        readlink_target = Path(os.readlink(local_sslg4))
+        require(readlink_target.name == source_cell["cell_id"]
+                and re.fullmatch(r"EJ(?:200|204|230)_x(?:m|p)\d+", readlink_target.name),
+                f"{source_cell['cell_id']}: invalid SSLG4 symlink target {readlink_target}")
```

Sources: `analysis/track_mechanism_20260915/prepare_campaign.py`, `git show f801d76 -- analysis/track_mechanism_20260915/prepare_campaign.py`, and the preserved `_BROKEN_20260916` directory.

The preserved timestamp for `EJ200_xm200/sslg4` is `2026-09-16 20:23:58.124548027 +0200`; commit `f801d76` is `2026-09-16T20:24:39+02:00`. The measured difference is approximately 41 seconds. This differs from the initial 15-second estimate; the filesystem timestamp and commit metadata are the recorded sources.

## 2. Regression test before the fix

The test was added at `analysis/track_mechanism_20260915/test_prepare_campaign_symlinks.py`. It prepares only a temporary campaign, then asserts that `readlink(cells/<cell>/sslg4).name` equals the receiving `cell_id`.

Against the uncorrected script it failed as follows:

```text
AssertionError: EJ200_xm200: unexpected target /home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/sslg4
```

This failure exercises the faulty source-selection logic rather than inspecting the contaminated campaign as a symptom.

## 3. Correction

`prepare_campaign.py` now defines material-specific validated sources:

| material | runtime source |
|---|---|
| EJ-200 | `/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/sslg4` |
| EJ-204 | `/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/sslg4` |
| EJ-230 | `/home/rrios/exec46_20260915/build_baseline/sslg4` |

For every cell, the script creates `runtime_by_cell/<cell_id>` as a symlink to the material-specific runtime, and `cells/<cell_id>/sslg4` as a symlink to that cell alias. Therefore the direct target path contains both material and position, while the resolved source remains the validated runtime for that material. The I2 `analysis_summary.json` gate remains mandatory: a corrected preparation requires `status == "PASS"`.

## 4. In-script verification gate

Immediately after creating each `cells/<cell_id>/sslg4` symlink, the script calls `os.readlink()` and requires both:

- the target basename equals the current `source_cell["cell_id"]`;
- the basename matches `EJ(?:200|204|230)_x(?:m|p)\d+`.

A mismatch raises `RuntimeError` with the cell and offensive target path before the loop can create another cell symlink. The positive test passed after this gate was added. Syntax compilation and `git diff --check` also passed.

## 5. Cleanup of the contaminated directory

The `_BROKEN_20260916` directory was not deleted. No process matching `detached_grid`, `ej200_bar_sim`, or the EXEC46 campaign path was running when checked with `pgrep -af`.

A `driver.lock` remains inside the preserved `_BROKEN_20260916` directory. The original path `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800` is absent because it was renamed to the evidence directory, so that preserved lock cannot interfere with preparation at the original path. The new v2 directory has no lock or simulation output.

## 6. v2 regeneration

A new campaign was prepared at:

`/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2`

The preparation used the baseline executable `/home/rrios/exec46_20260915/build_baseline/ej200_bar_sim` and the validated I2 corrected source argument. Preparation completed with 21 cells and `EXECUTED=0`. No ROOT from the contaminated campaign was reused.

The I2 gate was checked by the script against:

`/home/rrios/ej200/analysis/track_mechanism_20260915/i2_bc404_validation/analysis_summary.json`

which declares `status: PASS`.

The validated MPT hashes recorded and independently recomputed are:

| material | file | SHA-256 |
|---|---|---|
| EJ-200 / OPSC-100 | `rIndex.txt` | `15d1f8cf5a62effd9a0f2f9bd1edaeb5164b94c2035a11cfe6a878d51680f6ed` |
| EJ-200 / OPSC-100 | `absLength.txt` | `82191c023ed343d9b0f0c6de09ec3f1be12627e4eed6219d8e048529d7755d3a` |
| EJ-204 / OPSC-101 | `rIndex.txt` | `f1c77ee162c767cd23be608e0f8081174263352e40a0e1f57b8ca97ccef8636f` |
| EJ-204 / OPSC-101 | `absLength.txt` | `9473735344b5bd883a748df0223129b692edb659a7bd9c18f8d32be450749660` |

EJ-230 remains `UNCORRECTED_CONSTANT_RINDEX_NO_MEASURED_ANALOG` and resolves uniformly to the baseline runtime.

## 7. Independent 21-cell symlink verification

The following table is the external check using `os.readlink()` on every v2 cell. The resolved source was also checked: EJ-200 resolves to the F4 runtime, EJ-204 to the I2 runtime, and EJ-230 to baseline. Every row passed the material-specific MPT check where applicable.

| cell | `readlink(cells/*/sslg4)` | resolved source | check |
|---|---|---|---|
| EJ200_xm200 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ200_xm200` | `/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/sslg4` | PASS |
| EJ200_xm500 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ200_xm500` | `/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/sslg4` | PASS |
| EJ200_xm650 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ200_xm650` | `/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/sslg4` | PASS |
| EJ200_xp0 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ200_xp0` | `/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/sslg4` | PASS |
| EJ200_xp200 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ200_xp200` | `/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/sslg4` | PASS |
| EJ200_xp500 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ200_xp500` | `/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/sslg4` | PASS |
| EJ200_xp650 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ200_xp650` | `/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/sslg4` | PASS |
| EJ204_xm200 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ204_xm200` | `/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/sslg4` | PASS |
| EJ204_xm500 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ204_xm500` | `/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/sslg4` | PASS |
| EJ204_xm650 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ204_xm650` | `/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/sslg4` | PASS |
| EJ204_xp0 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ204_xp0` | `/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/sslg4` | PASS |
| EJ204_xp200 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ204_xp200` | `/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/sslg4` | PASS |
| EJ204_xp500 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ204_xp500` | `/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/sslg4` | PASS |
| EJ204_xp650 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ204_xp650` | `/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/sslg4` | PASS |
| EJ230_xm200 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ230_xm200` | `/home/rrios/exec46_20260915/build_baseline/sslg4` | PASS |
| EJ230_xm500 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ230_xm500` | `/home/rrios/exec46_20260915/build_baseline/sslg4` | PASS |
| EJ230_xm650 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ230_xm650` | `/home/rrios/exec46_20260915/build_baseline/sslg4` | PASS |
| EJ230_xp0 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ230_xp0` | `/home/rrios/exec46_20260915/build_baseline/sslg4` | PASS |
| EJ230_xp200 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ230_xp200` | `/home/rrios/exec46_20260915/build_baseline/sslg4` | PASS |
| EJ230_xp500 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ230_xp500` | `/home/rrios/exec46_20260915/build_baseline/sslg4` | PASS |
| EJ230_xp650 | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/runtime_by_cell/EJ230_xp650` | `/home/rrios/exec46_20260915/build_baseline/sslg4` | PASS |

All 21 rows passed the independent check. The v2 tree contains no ROOT output from a production attempt.

## 8. v2 dry-run

Command executed:

```bash
python3 analysis/sigma_t/orchestration/detached_grid.py dry-run \
  --directory /home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2
```

Recorded result:

| field | value |
|---|---:|
| status | `PASS` |
| job count | 21 |
| readable macros | 21 |
| fresh attempt output directories | 21 |
| diagnostics | `false` |
| timeout | `null` |
| execution marker | `EXECUTED=0` |

The independent pre-dry-run check found zero ROOT outputs. No `launch` command was executed.

## 9. Launch hold point

The exact launch command, printed but not executed, is:

```bash
python3 /home/rrios/ej200/analysis/sigma_t/orchestration/detached_grid.py launch \
  --directory /home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2
```

Awaiting explicit approval to launch v2.
