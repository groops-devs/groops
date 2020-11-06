#!/bin/bash
../bin/groops --doc .

echo "make Latex"
cd latex
pdflatex -interaction=nonstopmode documentation.tex >/dev/null
pdflatex -interaction=nonstopmode documentation.tex >/dev/null
pdflatex -interaction=nonstopmode documentation.tex | grep -a1 -e "^!" -e "Warning"
cd ..
ln -f -s latex/documentation.pdf

echo "make source code docu (doxygen)"
cd source
doxygen Doxyfile >/dev/null
cd ..
