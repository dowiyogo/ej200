# Reglas Permanentes de Operación para Agentes (`t0minidaq` / `dowiyogo/ej200`)

**Lectura obligatoria antes de ejecutar cualquier comando o editar cualquier archivo en `/home/rrios/`.**

Este host (`t0minidaq`) aloja el repositorio canónico `/home/rrios/ej200` (`dowiyogo/ej200`, simulación Geant4 11.x + SSLG4-OPSim y análisis ROOT/PyROOT del SHiP Timing Detector T0) junto con más de 1 TB de campañas de simulación y datos de tesis. Cualquier agente que opere en `/home/rrios/` está sujeto sin excepción a las siguientes reglas operativas.

---

## 1. Documentos de Referencia Obligatorios (Leer Antes de Actuar)

Antes de proponer, clasificar, mover o editar archivos, consulta siempre las fuentes de verdad ya auditadas en `/home/rrios/ej200`:

1. [`docs/reports/AUDITORIA_REPO_20260925.md`](file:///home/rrios/ej200/docs/reports/AUDITORIA_REPO_20260925.md) — Censo completo de árboles Git, discrepancias de material en campañas `.root`, genealogía de ramas y trazabilidad de figuras.
2. [`docs/execution_logs/REORG_20260925.md`](file:///home/rrios/ej200/docs/execution_logs/REORG_20260925.md) — Registro de la reorganización en 8 pasos, tabla de hashes SHA-256 de archivos `.root` resguardados, criterios formales de clasificación y bitácora de housekeeping.
3. [`docs/catalogo_macros.md`](file:///home/rrios/ej200/docs/catalogo_macros.md) — Catálogo de las 115 macros ROOT (`.C`), su estado de acoplamiento (`#include`, `.L`, invocación por `rebuild*.sh` o `.tex`) y los 14 módulos Python internos importados por otros scripts.
4. [`docs/reports/INVENTARIO_HUERFANOS_20260925.md`](file:///home/rrios/ej200/docs/reports/INVENTARIO_HUERFANOS_20260925.md) — Clasificación por contenido real de archivos sueltos, reportes, bitácoras, presentaciones, figuras sin sidecars completos y artefactos huérfanos.

---

## 2. Reglas Estrictas sobre Archivos, Rutas y Control de Versiones

1. **Verificación previa y posterior con `grep -rln` (prohibido `sed` masivo ciego):**
   - Antes de mover o renombrar cualquier archivo, ejecuta `grep -rln "<nombre_actual>"` sobre **todo** el árbol relevante (excluyendo binarios pesados como `*.root`, `*.npz`, `*.pdf`, `*.png`, `.git/`, `build*/`, `.venv/`).
   - Enumera cada referencia encontrada y edita cada archivo consumidor **uno por uno** con contexto explícito. Jamás ejecutes reemplazos masivos a ciegas con `sed -i` o equivalentes sobre múltiples archivos.
   - Tras aplicar el movimiento y actualizar las referencias, vuelve a correr exactamente el mismo `grep -rln` para confirmar **cero resultados** fuera de bitácoras históricas.

2. **Tag de respaldo previo y commits atómicos con trazabilidad `archivo:línea`:**
   - Antes de iniciar cualquier bloque de cambios en `/home/rrios/ej200`, crea un tag de respaldo sobre `HEAD` (p. ej., `git tag pre-<tarea>-<AAAAMMDD>`).
   - Después de cada cambio lógico autocontenido, realiza un **commit atómico** cuyo cuerpo detalle explícitamente cada modificación en formato `archivo:línea (antes -> después)`.

3. **Prohibición absoluta de tocar archivos `.root` sin autorización explícita:**
   - **NUNCA** muevas, renombres, sobrescribas ni borres ningún archivo `.root` (ni de simulación cruda, ni derivado de análisis, ni sidecar de figura) sin aprobación explícita del usuario en el chat para ese archivo específico. Los `.root` no están versionados en Git y su pérdida o corrupción es irreversible.

4. **Operaciones reservadas exclusivamente al usuario (`git push`, `git worktree remove`, `rm -rf`):**
   - **NUNCA** ejecutes `git push`, `git worktree remove`, `git branch -D` ni `rm -rf` sobre clones, worktrees o directorios de campaña.
   - Cuando una tarea requiera alguna de estas operaciones, imprime el bloque exacto de comandos verificados y detente para que el usuario los ejecute manualmente desde su terminal.

5. **Protección de macros ROOT/CINT (`.C`) y acoplamiento nombre-función:**
   - **NUNCA** renombres una macro `.C` sin consultar antes [`docs/catalogo_macros.md`](file:///home/rrios/ej200/docs/catalogo_macros.md) y verificar con `grep -rn` si tiene acoplamiento crítico (otros archivos que la cargan vía `#include`, `gROOT->ProcessLine(".L ...")`, scripts `rebuild*.sh` o decks `.tex`).
   - En ROOT/CINT y ACLiC, el nombre del archivo `.C` debe coincidir exactamente con el nombre de la función de entrada declarada adentro (`void nombre_archivo(...)`). Renombrar el archivo sin actualizar simultáneamente la firma interna y todos sus invocadores rompe la macro de forma silenciosa.

6. **Prohibido adivinar o inventar nombres de archivos o rutas:**
   - **NUNCA** escribas ni uses en un comando una ruta o nombre de archivo que no hayas verificado previamente con una herramienta de lectura real (`find`, `ls`, `stat`, `git status`). Si no tienes certeza absoluta del nombre exacto en disco, confírmalo primero con `find` o `ls`.

7. **Preservación estricta de `mtime` en reubicaciones:**
   - Para mover archivos dentro o hacia el repositorio, usa siempre `git mv` o `mv` (nunca `cp` seguido de `rm` sin preservar metadatos).
   - El timestamp de modificación (`mtime`) es parte de la evidencia cronológica de las campañas `EXEC_N`; un `mtime` nuevo solo es legítimo cuando el contenido del archivo fue editado funcionalmente (por ejemplo, al actualizar una ruta o un `import`).

---

## 3. Reglas de Rigor Científico y Clasificación

8. **Identificación del material real de archivos `.root` por evidencia directa:**
   - **NUNCA** infieras el material centellador de un archivo `.root` únicamente por el nombre del directorio o del archivo (existen precedentes documentados de carpetas llamadas `ej200_*` o `ej230_*` cuyos `.root` fueron corridos con otro material).
   - Determina siempre el material real cruzando evidencia directa: el `run.mac` / `campaign.json` de la corrida (`/det/scintillator` o `/det/scintillatorMaterial`), el `git log` del commit del binario que lo generó y el contenido de los `TTree`s.
   - **Atención física crítica:** `EJ-228` (`OPSC-105`) y `EJ-230` (`OPSC-106`) comparten exactamente el mismo espectro de emisión con pico en **$391\text{ nm}$** según las hojas de datos de Eljen Technology; por lo tanto, el branch `wl_nm` por sí solo **NO** distingue `EJ-228` de `EJ-230` (debe verificarse el código OPSC, el tiempo de decaimiento $\tau_d$ o la geometría barra vs cilindro).

9. **Regla de no-adivinar ante casos ambiguos:**
   - Si la clasificación, el origen o el impacto de tocar un archivo es genuinamente ambiguo (no hay evidencia concluyente ni en el contenido del archivo ni en el historial Git), **detente en ese archivo**, no adivines ni apliques reglas por analogía, regístralo en la tabla de **"casos dudosos"** de la bitácora y continúa con el resto de la tarea.

10. **Criterios formales para clasificar código y campañas ("activo" vs "histórico/obsoleto"):**
    - Antes de etiquetar cualquier script, macro o directorio de campaña como `"activo"`, `"vivo"`, `"histórico"` o `"obsoleto"`, aplica estrictamente las definiciones y reglas formales ya establecidas en [`docs/execution_logs/REORG_20260925.md`](file:///home/rrios/ej200/docs/execution_logs/REORG_20260925.md) y [`docs/reports/AUDITORIA_REPO_20260925.md`](file:///home/rrios/ej200/docs/reports/AUDITORIA_REPO_20260925.md). No inventes criterios ad-hoc en cada sesión.

11. **Alcance completo de `/home/rrios/` y verificación obligatoria de carpetas `build*`:**
    - Cualquier auditoría futura debe barrer TODO `/home/rrios/` con `find -maxdepth 1`, no solo los clones Git (ver [`docs/reports/AUDITORIA_CAPA_SUELTA_20260925.md`](file:///home/rrios/ej200/docs/reports/AUDITORIA_CAPA_SUELTA_20260925.md)). Nunca asumir que una carpeta `build*` es 100% regenerable sin antes buscar `.root`/`.csv`/`.json` adentro y verificar su contenido real -- ya hubo campañas reales y hasta una carpeta build activa (`build_baseline`) escondidas ahí. Antes de borrar cualquier `build*`, confirmar con `grep` que ningún script activo la referencia por ruta.

---

> **Nota de referencia cruzada:** Dentro del repositorio Git existe una versión breve en [`/home/rrios/ej200/AGENTS.md`](file:///home/rrios/ej200/AGENTS.md) que reenvía a este documento completo e indexa los reportes clave de `docs/`.
