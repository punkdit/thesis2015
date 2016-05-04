

out: thesis
	open thesis.pdf 

thesis: thesis.tex refs.bib
	pdflatex thesis.tex
	bibtex thesis
	pdflatex thesis.tex
	pdflatex thesis.tex


supplement.tex: ../fibonacci/supplement.tex 
	./build.py ../fibonacci supplement.tex

repr.tex: ../spectrum2015/repr.tex 
	./build.py ../spectrum2015 repr.tex




spectrum: spectrum.tex refs.bib
	pdflatex spectrum.tex
	bibtex spectrum
	pdflatex spectrum.tex
	pdflatex spectrum.tex
	open spectrum.pdf 


