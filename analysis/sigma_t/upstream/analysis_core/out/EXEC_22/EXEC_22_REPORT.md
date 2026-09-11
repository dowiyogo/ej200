# EXEC_22_REPORT — 2026-06-21

## 1. Veredicto: COMPLETADO-CON-FLAGS (fix parcial)

El fix TIR (dielectric_dielectric) fue implementado y auditado. La auditoria de velocidad
PASA. Sin embargo, las 3 condiciones de T3 no se cumplen completamente en el centro debido
a que la recuperacion de fotones no-TIR mediante REFLECTIVITY en el MPT es insuficiente en
Geant4 11.04. Se documenta el problema y la ruta al fix completo.

---

## 2. Velocidad (T2): PASA ✓

- Fraccion superluminica: 0.0000% (0 fotones con v > c/n en ninguna posicion)
- Max v/(c/n): 0.9721 (bien por debajo de 1.0)
- Sin gap de aire → sin fotones en aire → sin riesgo superluminico
- T4 EndTop mini-sim: tambien PASA (max v/(c/n) = 0.936)

---

## 3. Fix aplicado (T1)

### ej200_endonly (exec21-optfix branch)
Commits: af5ddb7 (exec21: dielectric_metal→dielectric_dielectric)
         9b8361f (exec22: REFLECTIVITY=0.95 para fotones no-TIR)

### ej200 EndTop build (exec22-endtop-optfix branch)
Commit: 76b582c (port del mismo fix a CreateBarSkinReflector)

### Cambios en Materials.cc:
```
surf->SetType(dielectric_dielectric);  // antes: dielectric_metal
surf->SetFinish(polished);
mpt->AddProperty("REFLECTIVITY", {1.5eV,6.5eV}, {0.95,0.95});
```

### Para re-simular la campana EndTop completa (NO ejecutar sin aprobacion de Rene):
```bash
cd /home/reriosto/SHiP/ej200
git checkout exec22-endtop-optfix
# Compilar y lanzar scan EndTop 5000 ev/pos:
# [analogo a t0minidaq_endtop_scan_20260618_204959 pero con el build corregido]
```

---

## 4. Condiciones T3 (3 pos, 500 ev)

| Condicion | Resultado | Comentario |
|----------|----------|-----------|
| (a) Npe realista | FALLA (centro: 0.4 PE, real: ~40 PE) | REFLECTIVITY no recupera no-TIR en Geant4 11.04 |
| (b) σ_END < 250 ps | FALLA (centro: ~985 ps) | Mismo problema |
| (c) Bimodal time structure | PASA (48% prompt @ x=-690) | TIR SI funciona! |

**Mejora real del fix**: perfil Npe MAS SUAVE cerca extremos. x=-600: 26.1 PE (antes ~0).
La TIR esta activa (condicion c pasa). El no-TIR sigue escapando sin recuperacion.

---

## 5. Por que falla el centro (diagnostico del fix parcial)

En Geant4 11.04, `REFLECTIVITY` en el MPT de una superficie `dielectric_dielectric`
unificada no overrides el comportamiento de Fresnel para fotones no-TIR de la manera
esperada. Los fotones con angulo < theta_c = 39.3° siguen transmitiendo al world volume
en lugar de reflejarse con probabilidad R=0.95.

**Consecuencia**: 23% de fotones por rebote siguen escapando.
Con 164 rebotes promedio (centro→END): (0.77)^164 ≈ 0 survival.

**Fix completo requerido**: geometria explicita:
1. Volumen de aire (gap ~0.1mm) rodeando la barra
2. Panel Mylar/ESR como volumen hijo del world con superficie `dielectric_metal` R=0.95-0.98
3. Fotones no-TIR que pasan a aire se reflejan en el Mylar y vuelven al centellador

---

## 6. T5 - Baseline END-only (scan 500 ev/pos, R=0.90)

Scan completado: 31 posiciones, 500 ev/pos, fix parcial (R=0.90 por defecto del script).
Perfil Npe: suave cerca extremos, colapso en centro. No fisico todavia.

---

## 7. T6 - TOP contribution (provisional/invalido)

σ_END no fisico → "aporte TOP 85%" de EXEC_20 sigue invalido.
Recalcular cuando fix completo este implementado.

GLS fix (EXEC_21 D6): σ_EndTop ≈ σ_TOP ≈ 76 ps (inv-varianza; unidades corregidas).

---

## 8. Commits y tags

- analysis_core master: EXEC_22-pre-t1 (tag), EXEC_22-final (tag)
- ej200_endonly exec21-optfix: af5ddb7, 9b8361f
- ej200 exec22-endtop-optfix: 76b582c

git push origin master  # (NO ejecutado)
