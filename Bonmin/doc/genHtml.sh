#!/bin/bash
for f in $@ ; 
do
NAME=`basename $f .tex`;
rm -rf $NAME;
awk -v NAME=$NAME '{sub(/FILENAME/,NAME);print}' Head.tex > tmp.tex;
latex tmp.tex;
latex2html -split 1 -no_navigation -info 0  tmp.tex;
awk '{sub(/tmp.css/,"../bonmin.css") ; sub(/>tmp</,">BONMIN Users Manual<") ; print}' tmp/index.html > toto.html;
mv toto.html tmp/index.html;
mv tmp $NAME;
#rm tmp.*;
done
