#!/bin/bash
curdir=$(pwd)
if [ $# != 4 ]
then
    echo "usage: $0 author date title notebook"
    exit -1
else
    author=$1
    date=$2
    title=$3
    notebook=$4
fi
mkdir -p tmp
cd tmp
echo """
((*- extends 'article.tplx' -*))

((* set cell_style = 'style_python.tplx' *))

((* block output_prompt *))
((* endblock output_prompt *))

((* block docclass *))
\documentclass{beamer}
((* endblock docclass *))

((* block margins *))
((* endblock margins *))

((* block title *))
    \title{$title}
((* endblock title *))

((* block date *))
    \date{$date}
((* endblock date *))

((* block author *))
    \author{$author}
((* endblock author *))

((* block packages *))
    ((( super() )))
((* endblock packages *))

% Use reveal.js tags to handle skipping/ new slides etc
((* block any_cell *))
    ((* if cell.metadata.get('slideshow',{}).get('slide_type','-').lower() in ['skip'] *))
    ((* elif cell.metadata.get('slideshow',{}).get('slide_type','-').lower() in ['slide'] *))
        \newpage
        \section{ ((( cell.source ))) }
        \frame{\sectionpage}
    ((* elif cell.metadata.get('slideshow',{}).get('slide_type','-').lower() in ['subslide'] *))
        \newpage
        ((( super() )))
    ((* elif cell.metadata.get('slideshow',{}).get('slide_type','-').lower() in ['fragment'] *))
        ((( super() )))
    ((* else *))
        ((( super() )))
    ((* endif *)) 
((* endblock any_cell *))
""" > notebook.tplx
ipython nbconvert --quiet --to=latex --post=PDF\
    --template="notebook.tplx" \
    --post=PDF "$curdir/$notebook"
# vim:ts=4:sw=4:expandtab
