#!/bin/bash

categories="bonmin filter ipopt"

mkdir -p html

for f in $categories ; do
  awk -v NAME=$f '{sub(/TEMPLATE/,NAME) ; print}' options_list_TEMPLATE.tex > options_list_$f.tex
  ./genHtml.sh options_list_$f.tex
  mv tmp.html options_list_$f.html
done

cp ../bonmin.css html/
