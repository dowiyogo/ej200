.PHONY: analyzer analyzer-run analyzer-clean exec12t-analysis exec12t-report exec12t-beamer exec12t-all exec12tb-figures exec12tb-tables exec12tb-beamer exec12tb-qa exec12tb-all

analyzer:
	$(MAKE) -C analysis

analyzer-run:
	$(MAKE) -C analysis run

analyzer-clean:
	$(MAKE) -C analysis clean

EXEC12T_OUT ?= $(shell ls -dt results/exec12t_* 2>/dev/null | head -1)
exec12t-analysis:
	python3 analysis/exec12t_timing_threshold_analysis.py analyze --output-dir "$(EXEC12T_OUT)"
exec12t-report:
	python3 analysis/exec12t_make_products.py --out "$(EXEC12T_OUT)"
exec12t-beamer:
	bash scripts/build_exec12t.sh "$(EXEC12T_OUT)"
exec12t-all: exec12t-analysis exec12t-beamer

EXEC12TB_OUT ?= $(shell ls -dt results/exec12tb_* 2>/dev/null | head -1)
exec12tb-figures:
	python3 analysis/exec12tb_figures.py "$(EXEC12TB_OUT)"
exec12tb-tables:
	python3 analysis/exec12tb_tables.py  "$(EXEC12TB_OUT)"
exec12tb-beamer:
	cd "$(EXEC12TB_OUT)/beamer" && \
	  lualatex -interaction=nonstopmode exec12tb_beamer.tex && \
	  lualatex -interaction=nonstopmode exec12tb_beamer.tex
exec12tb-qa:
	@pdfinfo "$(EXEC12TB_OUT)/beamer/exec12tb_beamer.pdf" | grep Pages
	@python3 -c "\
import hashlib,pathlib; \
figs=list(pathlib.Path('$(EXEC12TB_OUT)/figures').glob('*.pdf')); \
shas={}; dups=0; \
[shas.setdefault(hashlib.sha256(f.read_bytes()).hexdigest(), []).append(f.name) for f in figs]; \
dups=sum(len(v)>1 for v in shas.values()); \
print(f'{len(figs)} figures, {dups} SHA duplicates'); \
exit(dups)"
exec12tb-all: exec12tb-figures exec12tb-tables exec12tb-beamer exec12tb-qa
