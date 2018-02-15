doc=main.tex
pdf=$(doc:.tex=.pdf)

all: $(pdf)

$(pdf): $(doc)
	pdflatex $<
	bibtex `basename -s .tex $<`
	pdflatex $<
	pdflatex $<

clean:
	rm -f *~ *.toc *.aux *.out *.log *.nav *.snm *.xml *-blx.bib *.bbl *.blg

proper: clean
	rm $(pdf)
