all:
	bibtex paper
	pdflatex paper
	make clean

clean:
	- rm paper.log
	- rm paper.aux
	- rm paper.bbl
	- rm paper.dvi
	- rm paper.blg
	- rm paper.out
	- rm -f paperfigures/*converted-to.pdf