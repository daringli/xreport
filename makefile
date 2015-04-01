tex = report2.tex

love: report2.pdf

anyways: 
	pdflatex report2.tex
	bibtex report2.aux
	makeindex report2.nlo -s nomencl.ist -o report2.nls
	pdflatex reportSimul2013.tex
	pdflatex reportSimul2013.tex

report2.pdf:  bibli.bib joelspropack.bst report2.tex titlepageabstract/titlepageabstract.tex preamble/preamble.tex titlepageabstract/acknowledge/acknowledgements.tex titlepageabstract/abstract/sammandrag.tex titlepageabstract/abstract/abstract.tex intro/intro.tex  
	pdflatex report2.tex
	bibtex report2.aux
	makeindex report2.nlo -s nomencl.ist -o report2.nls
	pdflatex report2.tex
	pdflatex report2.tex

war: temp
	rm report2.pdf
  
temp:
	rm report2.aux
	rm report2.bbl


#pdflatex nomenclature-test.tex
#pdflatex nomenclature-test.tex
#pdflatex nomenclature-test.tex
