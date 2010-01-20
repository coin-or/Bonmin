OUTPUTS = Install.done
OUTPUTS += Intro.done
OUTPUTS += Obtain.done
OUTPUTS += bib.done
OUTPUTS += options_list.done
OUTPUTS += options_set.done
OUTPUTS += use.done

all: $(OUTPUTS) Head.tex

%.done : %.tex Head.tex
	./genHtml.sh $<
	touch $@

install_www:
	cp -rv bib Install Intro Obtain options_list options_set use ../../../../html/
