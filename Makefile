
.PHONY: writeup view clean

writeup: writeup.tex
	latexmk -xelatex -use-make writeup.tex
	open writeup.pdf

writeup.tex: writeup.md
	pandoc -s -S  --template=latex.template --number-sections -o writeup.tex writeup.md

view:
	open writeup.pdf

clean:
	latexmk -CA