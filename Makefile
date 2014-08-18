all: open-science.html open-science.pdf

clean:
	rm -f open-science.html open-science.pdf

install-deps:
	brew install pandoc

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

# Dependencies

open-science.pdf: homebrew-bioinformatics.png

# Rules

%.html: %.md
	pandoc -st revealjs -V theme:sky $< \
		|sed 's/history: true/slideNumber: true/' >$@

%.pdf: %.md
	pandoc -t beamer -o $@ $<

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
