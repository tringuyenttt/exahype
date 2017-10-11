.IGNORE:

# keep a list of texfiles here
TARGETS = exahype-advanced-hydro

# call this makefile like:
#   make exahype-z4         # to make an indivdual PDF
#   make all                # to make all PDFs (default)

pdflatex = pdflatex -interaction=nonstopmode

all: $(TARGETS)

$(TARGETS):
	$(pdflatex)     "$@"
	#bibtex         "$(TARGET)"
	$(pdflatex)     "$@"
	$(pdflatex)     "$@"

.PHONY: clean vclean

clean:
	/bin/rm -f *.dvi *.log *.aux *.toc *.bbl *blg *.end *.out *.4ct *.4tc *.idv *.lg *.tmp *.xref *.synctex.gz

vclean: clean
	/bin/rm -f $(TARGETS:=.pdf)  *~
	/bin/rm -rf build/

