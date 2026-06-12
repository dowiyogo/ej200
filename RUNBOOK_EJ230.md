# RUNBOOK — EXEC_13 / EJ-230 campaign

Branch: `feat/ej230-sslg4`  
Scintillator: EJ-230 (SSLG4 OPSC-106), 9700 ph/MeV, τ_r=0.5 ns, τ_d=1.5 ns, att.=120 cm  
Geometry: BYTE-IDENTICAL to feat/endtop-sslg4 (16 End + 70 Top SiPMs, same seeds, same N=2000 ev/pos)

---

## En t0minidaq (usuario: rrios)

```bash
# 1. Clonar y configurar
cd /home/rrios
git clone https://github.com/dowiyogo/ej200 ej230
cd ej230
git checkout feat/ej230-sslg4

# 2. Compilar (GDML OFF porque t0minidaq no tiene soporte Xerces-C)
mkdir build && cd build
cmake .. -DENABLE_GDML_TEST=OFF
make -j24
ctest   # todos los tests excepto GDML deben pasar

# 3. Test al centro (validación rápida, 200 eventos)
../scripts/run_center_t0minidaq_24t.sh

# 4. Scan largo (31 posiciones × 2000 eventos)
../scripts/run_scan_t0minidaq_24t.sh

# 5. Análisis completo (ROOT → CSV → figuras → Beamer .tex)
../scripts/run_analysis_t0minidaq.sh

# 6. Descargar resultados en MSI
# El script imprime el tarball y el scp exacto al final, p. ej.:
#   scp rrios@t0minidaq.cern.ch:/home/rrios/ej230/results_ej230_YYYYMMDD_HHMM.tar.gz .
```

**Notas t0minidaq:**
- `ENABLE_GDML_TEST=OFF` es obligatorio (GDML no disponible en esta máquina).
- Las dependencias Python (`uproot`, `numpy`, `pandas`, `matplotlib`, `scipy`, `awkward`)
  se instalan automáticamente si faltan (pip --user dentro del script de análisis).
- Todo el análisis corre headless (`matplotlib.use("Agg")`; ROOT no se invoca).
- Los resultados se escriben bajo `$HOME/ej230/results_ej230/` (configurable con `RESULTS_DIR`).
- Nunca se sobreescriben ROOT existentes: cada reanudación del scan crea un subdirectorio
  con timestamp bajo `data/`.

---

## En MSI (usuario: reriosto)

```bash
# Build estándar (GDML ON por defecto — Geant4 11.4.0 + Xerces-C presentes)
cd /home/reriosto/SHiP/ej230
mkdir build && cd build
cmake ..          # ENABLE_GDML_TEST=ON por defecto
make -j16
ctest             # todos los tests deben pasar, incluido GDML

# Ejecución local (si se desea validar sin t0minidaq)
../scripts/run_center_msi_16t.sh
../scripts/run_scan_msi_16t.sh      # ← scan largo; 16 threads; usa $HOME/ej230/results_ej230/

# Compilar Beamer desde árbol descargado de t0minidaq
tar -xzf results_ej230_YYYYMMDD_HHMM.tar.gz
RESULTS_DIR=$PWD/results_ej230 ../scripts/build_report_msi.sh
# → report/exec13_ej230_report_full.pdf
```

**Notas MSI:**
- Geant4 11.4.0 + GDML disponibles; todos los CTests activos.
- `sslg4_properties_check` verifica yield=9700/MeV, rise=0.5 ns, decay=1.5 ns, att.=120 cm.
- Para regenerar solo el Beamer (sin nueva simulación):
  ```bash
  python3 scripts/make_beamer_ej230.py
  cd analysis/exec07 && pdflatex -interaction=nonstopmode exec13_ej230_report_full.tex
  ```

---

## Variables de entorno

| Variable | Default | Descripción |
|----------|---------|-------------|
| `RESULTS_DIR` | `$HOME/ej230/results_ej230` | Directorio raíz de resultados |
| `BUILD_DIR` | `$REPO/build` | Directorio de compilación |
| `EVENTS_PER_POINT` | `2000` | Eventos por posición en el scan |
| `EVENTS_CENTER` | `200` | Eventos en el run de validación central |

---

## Estructura de resultados (`RESULTS_DIR/`)

```
results_ej230/
├── data/                    # ROOT files (photon_hits_x{x}mm.root)
│   └── session_YYYYMMDD_*/  # Subdirectorios por sesión (nunca sobreescribir)
├── csv/                     # exec13_tN_summary.csv, ...
├── figures/                 # exec13_tN_{x}mm.png, exec13_tN_summary.png, ...
├── tables/                  # exec13_tN_t4_table.tex, exec13_macros.tex
├── report/                  # exec13_ej230_report_full.tex, .pdf
└── logs/                    # run_x{x}mm.log, exec13_tN_analysis.log, ...
```

---

## Verificaciones

| Verificación | Cuándo | Cómo |
|-------------|--------|------|
| Material properties | `ctest` | `sslg4_properties_check` (yield, rise, decay, att.) |
| Balance fotones | `ctest` | `endtop_balance_smoke` |
| Configuración readout | `ctest` | `readout_config_check` |
| GDML (solo MSI) | `ctest` | `export_endtop_gdml` + `endtop_gdml_check` |
| Paridad Beamer | `make_beamer_ej230.py` | frame count y títulos vs exec12 (119 frames) |

---

## Discrepancias documentadas

Ver [`docs/EJ230_material_audit.md`](docs/EJ230_material_audit.md) para:
- Tabla datasheet vs SSLG4 vs Geant4 para todos los parámetros EJ-230.
- `HOOK_PDE_SPECTRAL`: la curva PDE(λ) del SiPM produce ~2% menos eficiencia
  efectiva para EJ-230 (391 nm) que para EJ-204 (408 nm). Documentado, sin corrección.
