Scan preliminar longitudinal EJ-204 - 2026-06-08

Configuracion:
- Geant4 11.4.0, rama fix/physics-baseline
- Muon negativo, 1 GeV, incidencia vertical
- EJ-204 activo; rise time finito habilitado
- Posiciones x [mm]: -650, -350, -100, 0, 100, 350, 650
- 200 eventos por posicion
- /sipm/jitterSigma 0 ns para estimar sigma_t intrinseco

Estimador temporal:
- End: (FPT_left + FPT_right)/2, requiriendo hit en ambos extremos
- Top: primer fotoelectron entre IDs 16-35
- Ajuste gaussiano de la distribucion por posicion; fallback RMS

Eficiencia:
- 100 * fotones detectados / fotones de centelleo generados
- Contadores master de RunAction al final de cada beamOn

Advertencia:
- sigma_t es intrinseco y no incluye SPTR ni jitter de electronica.
