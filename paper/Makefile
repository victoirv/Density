COMP=pdflatex
BIB=bibtex
TEX=paper

all: 
	$(COMP) $(TEX).tex
	$(BIB) $(TEX)
	$(COMP) $(TEX).tex
	$(COMP) $(TEX).tex

clean:
	rm -f $(TEX).synctex.gz
	rm -f $(TEX).aux 
	rm -f $(TEX).bbl
	rm -f $(TEX).blg
	rm -f $(TEX).log
	rm -f $(TEX).out
	rm -f $(TEX).{ps,pdf,log,aux,out,dvi,bbl,blg,lof,lot,tks,toc}
	rm -f *.{log,aux} #Avoid removing ps and pdf, but get rid of chapter build files

diff:
	#checkout f26fb4b as that's pre-submission, compare to newest
	git checkout f26fb4b paper.tex
	mv paper.tex paper.old
	git checkout paper.tex
	latexdiff paper.old paper.tex > paper.diff.tex
	$(COMP) paper.diff.tex
	$(BIB) paper.diff
	$(COMP) paper.diff.tex
	$(COMP) paper.diff.tex
