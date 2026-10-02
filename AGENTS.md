# AGENTS.md

Este repositorio forma parte del ecosistema SHiP T0. **Las reglas obligatorias están en:**

- MSI: `/home/reriosto/SHiP/AGENTS.md`
- t0minidaq: `/home/rrios/AGENTS.md`

Léelas completas antes de tocar cualquier archivo. Las más importantes:

1. No crear clones ni copias de este repo. Para trabajar en paralelo: `git worktree add /tmp/<tarea>` y eliminarlo al terminar.
2. No crear archivos ni carpetas sueltos en la raíz del repo.
3. No cambiar fechas (`cp -p`, `rsync -a`; nunca `git archive` para transportar datos).
4. Para verificar o regenerar, trabajar en una copia en `/tmp`. Nunca crear symlinks a archivos.
5. No borrar con `rm -rf`. No dejar carpetas vacías.
6. Toda figura nueva lleva `.meta.json` y `.csv` en la misma carpeta.
7. Los datos (`.root`) no van a git. Viven en `t0minidaq:/home/rrios/ej200/data/ej200_campaigns/`, con `PROVENANCE.md`.
8. Nunca abrir, imprimir ni subir `.vscode/mcp.json` ni archivos con claves.
9. Nada de bucles `while kill -0` ni de subagentes como verificadores. Reportar salidas literales.

<!-- Reglas específicas de este repo (si las hay) van debajo de esta línea. -->

## Reglas específicas de este repo (texto anterior, sin cambios)

# Reglas de Operación para Agentes (`dowiyogo/ej200`)

Este repositorio se rige por las reglas permanentes de operación para agentes definidas en [`/home/rrios/AGENTS.md`](file:///home/rrios/AGENTS.md) (un nivel arriba de esta carpeta). **Léelas completas ANTES de mover, renombrar, borrar, editar o ejecutar cualquier cosa aquí dentro.**

## Referencias específicas de este repositorio (consultar según la tarea)

- [`docs/execution_logs/REORG_20260925.md`](file:///home/rrios/ej200/docs/execution_logs/REORG_20260925.md) — Historial completo de la reorganización (qué archivos se movieron, qué decisiones se tomaron, hashes SHA-256 de `.root` y qué quedó pendiente).
- [`docs/catalogo_macros.md`](file:///home/rrios/ej200/docs/catalogo_macros.md) — Catálogo de macros ROOT (`.C`) y módulos Python internos, indicando explícitamente qué macros NO se pueden renombrar por acoplamiento.
- [`docs/reports/AUDITORIA_REPO_20260925.md`](file:///home/rrios/ej200/docs/reports/AUDITORIA_REPO_20260925.md) — Censo técnico original del repositorio, discrepancias de material y trazabilidad de campañas y figuras.
- [`docs/reports/INVENTARIO_HUERFANOS_20260925.md`](file:///home/rrios/ej200/docs/reports/INVENTARIO_HUERFANOS_20260925.md) — Clasificación detallada de archivos sueltos y huérfanos por tipo de contenido real.
