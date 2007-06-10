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

