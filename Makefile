RM ?= $(shell dirname `which rm`)
PREFIX ?= $(shell dirname `which pdflatex`)

LATEX = perl texfot.pl $(PREFIX)/pdflatex -interaction=nonstopmode -halt-on-error -file-line-error
BIBTEX = $(PREFIX)/bibtex
DETEX = $(PREFIX)/detex

LINT = $(PREFIX)/chktex
LINT_OPTIONS = -q

.SUFFIXES: .tex .dvi .eps .ps .pdf

.PHONEY: all clean lint wc

MAIN = main
FIGDIR = figures
CHAPDIR = Appendix Communication gfc Lander Orbiter Penetrator Theory

FILES = $(notdir $(wildcard *.tex)) $(foreach dir,$(CHAPDIR),$(dir)/*.tex) bibliography/biblio.bib
FIGURES = $(shell find $(FIGDIR) -iname '*.png') $(shell find $(FIGDIR) -iname '*.jpg') $(shell find $(FIGDIR) -iname '*.pdf') $(shell find $(FIGDIR) -iname '*.eps') $(shell find $(FIGDIR) -iname '*.gif')

all: $(MAIN).pdf

$(MAIN).pdf: $(FILES) $(FIGURES)
	$(LATEX) $*.tex;
	$(BIBTEX) $*;
	$(LATEX) $*.tex 1>/dev/null;
	$(LATEX) $*.tex;

lint:
	$(LINT) $(LINT_OPTIONS) $(FILES) 2>/dev/null

clean:
	$(RM) -f *.aux \
		$(foreach dir,$(CHAPDIR),$(dir)/*.aux) \
		$(MAIN).log $(MAIN).dvi $(MAIN).ps $(MAIN).blg $(MAIN).bbl \
		$(MAIN).lot $(MAIN).lol $(MAIN).lof $(MAIN).toc $(MAIN).tdo $(MAIN).out $(MAIN).pdf

# Count words in the book
wc:
	@echo
	@echo "Current word count: "
	@$(DETEX) $(MAIN).tex | wc -w | sed 's/^[ \t]*/    /'
	@echo "Current page count: "
	@$(DETEX) $(MAIN).tex | tr -d ' \t\r\n' | echo "scale=3; `wc -m`/2200" | bc -l | sed 's/^[ \t]*/    /'
	@echo
