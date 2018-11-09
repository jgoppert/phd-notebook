#!/usr/bin/env bash
file=$1
filepath="${PWD}/${file}"
filename=$(basename "${file}")
filename="${filename%.*}"

mkdir -p tmp
cd tmp
echo """((*- extends 'article.tplx' -*))
((*- set cell_style = 'style_python.tplx' -*))
((*- block input -*))
    ((* if cell.metadata.suppress_input is defined and
    	cell.metadata.suppress_input == 1 *))
    ((* else *))
        (((super())))
    ((* endif *))
((*- endblock input -*))
((*- block output -*))
    ((* if cell.metadata.suppress_output is defined and
    	cell.metadata.suppress_output == 1 *))
    ((* else *))
        (((super())))
    ((* endif *))
((*- endblock output -*))

""" > my_latex.tplx
ipython nbconvert --quiet --to latex --post PDF \
	--template my_latex "${filepath}" && \
evince "${filename}.pdf" &
