all: open-science.html

clean:
	rm -f open-science.html

install-deps:
	brew install pandoc

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

# Dependencies

open-science.html: homebrew-bioinformatics.png

# Rules

%.html: %.md %.revealjs
	pandoc -st revealjs --template $* -V theme:sky -o $@ $<

# Render a graph of the dependencies of bioinformatics software in Homebrew-science

bioinformatics.tsv: /usr/local/Library/Taps/homebrew/homebrew-science/*
	grep -wl doi $^ |sed 's/.*\///;s/\.rb$$//' >$@

homebrew.gv:
	brew graph >$@ 2>/dev/null

homebrew-%.gv: homebrew.gv %.tsv
	(echo 'digraph $* {' && \
		grep -Ff $*.tsv $< && \
		echo '}') >$@

%.png: %.gv
	unflatten -fc9 -l9 $< |dot -Tpng >$@
