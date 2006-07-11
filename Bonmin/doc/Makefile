FILES = BOUM_Intro/index.html \
        BOUM_Obtain/index.html \
        BOUM_Install/index.html \
        BOUM_use/index.html \
        BOUM_options/index.html

all: $(FILES)
	@echo "html generated"

%/index.html: %.tex
	./genHtml.sh $<

clean: 
	rm -rf $(FILES)
