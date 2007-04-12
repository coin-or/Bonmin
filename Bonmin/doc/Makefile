OUTPUTS = BOUM_Install.done
OUTPUTS += BOUM_Intro.done
OUTPUTS += BOUM_Obtain.done
OUTPUTS += BOUM_bib.done
OUTPUTS += BOUM_options_list.done
OUTPUTS += BOUM_options_set.done
OUTPUTS += BOUM_use.done

all: $(OUTPUTS) Head.tex

%.done : %.tex Head.tex
	./genHtml.sh $<
	touch $@

