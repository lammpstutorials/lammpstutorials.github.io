#!/bin/bash
set -e

var1=$1
var2=$2

# convert the sphinx bibliography to ebook
if [[ "$var1" == "full" ]];
then
    python3 sphinx-to-ebook.py
fi
# compile pdf
pdflatex lammps-tutorials-ebook.tex
# if full keyword is given, recompile pdf
if [[ "$var1" == "full" ]];
then
    pdflatex lammps-tutorials-ebook.tex
fi

mv lammps-tutorials-ebook.pdf _lammps-tutorials-ebook.pdf 
pdftk logo/first-page.pdf _lammps-tutorials-ebook.pdf  cat output lammps-tutorials-ebook.pdf
rm _lammps-tutorials-ebook.pdf 

# recompile figures
# if [[ "$var2" == "fig" ]];
# then
#    basedir="$PWD"
    # collect all figure paths
    # find "../docs/inputs/level1/" -name 'plot_*.ipynb' > notebook-figures.txt
    #find "../docs/inputs/" -name 'plot_*.ipynb' > notebook-figures.txt
#    while read -r line
#    do
#        jupyter-nbconvert --to script $line
#        dir="$(dirname "${line}")"
#        ipynb="$(basename "${line}")"
#        name=`echo "$ipynb" | cut -d'.' -f1`
#        py=$name'.py'
#        cd $dir
#        python3 $py || echo "$line python failled"
#        cd $basedir
#    done < notebook-figures.txt
#fi