#!/bin/bash
mkdir -p html
for f in $@ ; 
do
NAME=`basename $f .tex`;
rm -rf $NAME;
awk -v NAME=$NAME '{sub(/FILENAME/,NAME);print}' Head.tex > tmp.tex;
latex tmp.tex;
latex tmp.tex;
latex tmp.tex;
tex4ht tmp.tex;
awk '{sub(/tmp.css/,"bonmin.css") ; sub(/>tmp</,">BONMIN Users Manual<") ; print}' tmp.html | sed -e's/##/#/g' > toto.html;
mv toto.html html/$NAME.html
#rm tmp.*
done
