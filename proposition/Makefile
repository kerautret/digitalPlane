doc=main.tex
pdf=$(doc:.tex=.pdf)
src=$(doc) section1.tex section2.tex section3.tex mybibliography.bib schedule.tikz tasks.tikz

all: $(pdf)

$(pdf): $(src)
	pdflatex $<
	bibtex `basename -s .tex $<`
	pdflatex $<
	pdflatex $<

clean:
	rm -f *~ *.toc *.aux *.out *.log *.nav *.snm *.xml *-blx.bib *.bbl *.blg

proper: clean
	rm $(pdf)
