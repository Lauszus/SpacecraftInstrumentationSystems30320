RM = /bin/rm

PREFIX ?= /usr/bin/
LATEX = $(PREFIX)pdflatex -halt-on-error -file-line-error
BIBTEX = $(PREFIX)bibtex
DETEX = $(PREFIX)detex

LINT = $(PREFIX)chktex
LINT_OPTIONS = -q

.SUFFIXES: .tex .dvi .eps .ps .pdf

.PHONEY: all clean lint wc

MAIN = main
FIGDIR = figures
#CHAPDIR = chapters

# Add your own .eps figures to this list.
#FIGURES = $(notdir $(wildcard $(FIGDIR)/*.png)) $(notdir $(wildcard $(FIGDIR)/*.jpg)) $(notdir $(wildcard $(FIGDIR)/*.pdf))

# Add your own LaTeX files to this list.
FILES = $(notdir $(wildcard *.tex)) bibliography/biblio.bib

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex $(FIGURES) $(FILES)
	$(LATEX) $*.tex;
	$(BIBTEX) $*;
	$(LATEX) $*.tex;
	$(LATEX) $*.tex;

lint:
	@ $(LINT) $(LINT_OPTIONS) *.tex $(CHAPDIR)/*.tex 2>/dev/null

clean:
	- $(RM) -f *.aux \
        $(CHAPDIR)/*.aux \
		$(MAIN).log $(MAIN).dvi $(MAIN).ps $(MAIN).blg $(MAIN).bbl \
		$(MAIN).lot $(MAIN).lol $(MAIN).lof $(MAIN).toc $(MAIN).pdf

# Count words in the book
wc:
	- @echo
	- @echo "Current word count: "
	- @$(DETEX) $(FILES) | wc -w
	- @echo
