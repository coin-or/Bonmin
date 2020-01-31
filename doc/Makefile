OUTPUTS = html/Install.html
OUTPUTS += html/Intro.html
OUTPUTS += html/Obtain.html
OUTPUTS += html/bib.html
OUTPUTS += html/options_list.html
OUTPUTS += html/options_set.html
OUTPUTS += html/use.html
OUTPUTS += html/bonmin.css

all: $(OUTPUTS) Head.tex

html/%.html : %.tex Head.tex bonmin.sty
	./genHtml.sh $<
	
html/bonmin.css : bonmin.css
	cp bonmin.css html/

install_www:
	cp -rv html/*.html ../../../../html/
