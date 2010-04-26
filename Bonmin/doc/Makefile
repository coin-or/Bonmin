OUTPUTS = html/Install.html
OUTPUTS += html/Intro.html
OUTPUTS += html/Obtain.html
OUTPUTS += html/bib.html
OUTPUTS += html/options_list.html
OUTPUTS += html/options_set.html
OUTPUTS += html/use.html

all: $(OUTPUTS) Head.tex

html/%.html : %.tex Head.tex bonmin.sty
	./genHtml.sh $<

install_www:
	cp -rv html/*.html ../../../../html/
