---
title: 'Open, reproducible science using Make, RMarkdown and Pandoc'
author: 'Shaun Jackman'
date: '2014-10-09'
---

Open, reproducible science
------------------------------------------------------------

### using Make, RMarkdown and Pandoc

Shaun Jackman [\@sjackman][]

2014-10-09 at [VanBUG][], Vancouver, Canada

[![Creative Commons Attribution License](images/cc-by.png)][cc-by]

[Fork me on GitHub!][]

[\@sjackman]: http://twitter.com/sjackman
[VanBUG]: http://www.vanbug.org/
[cc-by]: http://creativecommons.org/licenses/by/4.0/
[Fork me on GitHub!]: https://github.com/sjackman/open-science

Shaun Jackman
------------------------------------------------------------

| [Genome Sciences Centre][], BC Cancer Agency
| Vancouver, Canada
| [\@sjackman][]
| [github.com/sjackman][]
| [sjackman.ca][]

![](images/sjackman.jpg)

[Genome Sciences Centre]: http://bcgsc.ca
[github.com/sjackman]: https://github.com/sjackman
[sjackman.ca]: http://sjackman.ca

Open and reproducible science
------------------------------------------------------------

+ Open science
+ Repeatable science
	- by you
	- by others
+ Reproducible science

Open science
================================================================================

Open science
------------------------------------------------------------

+ Publish all research outputs
+ Archive manuscripts
+ Publish papers in open-access journals
+ Sign peer reviews
+ Participate in public discussion, like [Twitter][]

[Twitter]: http://twitter.com/

Publish all research outputs
------------------------------------------------------------

+ Slides on [SpeakerDeck][] or [slideshare][]
+ Posters on [figshare][]
+ Code and data on [GitHub][]
+ Data in a plain-text format, like TSV
+ Archive manuscript on [bioRxiv][] or [arXiv][]

[SpeakerDeck]: https://speakerdeck.com/
[slideshare]: http://www.slideshare.net/
[figshare]: http://figshare.com/
[GitHub]: https://github.com/
[bioRxiv]: http://biorxiv.org/
[arXiv]: http://arxiv.org/

Reproducible science
================================================================================

Repeatable science
------------------------------------------------------------

Given the same data and code&hellip;

. . .

+ Reproduce the same results
+ At least by yourself, this should be the minimum bar
+ Hopefully repeatable by others as well

Reproducible science
------------------------------------------------------------

Given the manuscript&hellip;

. . .

Another scientist can

+ Repeat the experiment
+ Analyse the data
+ Draw the same conclusion

Repeatable vs. reproducible science
------------------------------------------------------------

Reproducibility is fundamental to science

. . .

+ Often we don't even accomplish repeatable science
+ So let's start there

Repeatable science
================================================================================

Managing software
------------------------------------------------------------

+ Install software using [Homebrew][] or [Linuxbrew][]
+ [Homebrew-science][] has tons of science software
+ To report versions of all software used, need only write&hellip;

> We used Linuxbrew to install the required software from Homebrew-science version 2014-08.

[Homebrew]: http://brew.sh
[Linuxbrew]: http://brew.sh/linuxbrew/
[Homebrew-science]: http://brew.sh/homebrew-science/

Homebrew navigates dependency hell
------------------------------------------------------------

![Dependencies of bioinformatics tools in Homebrew](images/homebrew-bioinformatics.png)

Publish data
------------------------------------------------------------

> Best way to set back your competitors is to release your #data. That way they have to analyze their data & *your* data

| C. Titus Brown [\@ctitusbrown][]
| [BOSC 2014 keynote][]
| [A History of Bioinformatics (in the Year 2039)][]

[\@ctitusbrown]: https://twitter.com/ctitusbrown
[BOSC 2014 keynote]: http://video.open-bio.org/video/1/a-history-of-bioinformatics-in-the-year-2039
[A History of Bioinformatics (in the Year 2039)]: http://www.slideshare.net/c.titus.brown/2014-bosckeynote

Version control
------------------------------------------------------------

+ git/[GitHub][] for (almost) everything!
+ Maybe not big, raw data
+ For experimental design data
+ For results and summary statistics
+ Data in a plain-text format, like TSV
+ GitHub [renders TSV][] pretty!

![GitHub renders TSV pretty!](images/GitHub-tsv.png)

------------------------------------------------------------

### GitHub [renders TSV][] pretty!

![GitHub renders TSV pretty!](images/GitHub-tsv.png)

[renders TSV]: https://help.github.com/articles/rendering-csv-and-tsv-data

A reproducible manuscript
================================================================================

One Makefile
------------------------------------------------------------

+ Downloads the data
+ Runs the command-line programs
+ Performs the statistical analyses using [R][]
+ and Generates the TSV tables
+ Renders the figures using [ggplot2][]
+ Renders the supplementary material using [RMarkdown][]
+ Renders the manuscript using [Pandoc][]

[R]: http://www.rstudio.com/
[ggplot2]: http://ggplot2.org/
[RMarkdown]: http://rmarkdown.rstudio.com/
[Pandoc]: http://johnmacfarlane.net/pandoc/

Turns this
------------------------------------------------------------

![[UniqTag Markdown][]](images/UniqTag-md.png)

[UniqTag Markdown]: https://github.com/sjackman/uniqtag-paper/blob/master/UniqTag.md

Into this
------------------------------------------------------------

![[UniqTag PDF][]](images/UniqTag-pdf.png)

[UniqTag PDF]: http://biorxiv.org/content/early/2014/08/01/007583.full.pdf

Workflow
================================================================================

------------------------------------------------------------

[Plain Text, Papers, Pandoc][] by [Kieran Healy][]

![I promise this is less insane than it appears](images/workflow-rmd-md.png)

[Plain Text, Papers, Pandoc]: http://kieranhealy.org/blog/archives/2014/01/23/plain-text/
[Kieran Healy]: http://kieranhealy.org/

Make is beautiful
------------------------------------------------------------

| Tell Make how to create one type of file from another
| and which files you want to create.

. . .

| Make looks at which files you have
| and figures out how to create the files that you want.

Make example
------------------------------------------------------------

```makefile
%.bam: %.sam
	samtools view -Sb $< >$@

%.sort.bam: %.bam
	samtools sort $< $*.sort

%.bam.bai: %.bam
	samtools index $<
```

```bash
touch hello.sam
make hello.sort.bam.bai
```

```bash
samtools view -Sb hello.sam >hello.bam
samtools sort hello.bam hello.sort
samtools index hello.sort.bam
```

Markdown for the manuscript
------------------------------------------------------------

Markdown is a plain-text typesetting language

```markdown
A header
========

A list:

+ This text is *italic*
+ This text is **bold**
```

### A header

A list:

+ This text is *italic*
+ This text is **bold**

RMarkdown
------------------------------------------------------------

+ [RMarkdown][] interleaves text with code in [R][]
+ Code that calculates summary statistics
+ Code that generates tables
+ Code that renders figures using [ggplot2][]
+ [RMarkdown][] is ideal for supplementary material

Pandoc
------------------------------------------------------------

| [Pandoc][] renders attractive documents and slides
| from plain-text typesetting formats

It converts between every format known (just about)

+ Markdown
+ HTML
+ LaTeX
+ PDF
+ ODT and docx (yes, really)

fin
================================================================================

Links
------------------------------------------------------------

| [SpeakerDeck][] | [slideshare][] | [figshare][]
| [GitHub][] | [bioRxiv][] | [arXiv][]
| [Homebrew][] | [Linuxbrew][] | [Homebrew-science][]
| [R][] | [ggplot2][] | [RMarkdown][] | [Pandoc][]
| [Plain Text, Papers, Pandoc][]

[\@sjackman][] | [github.com/sjackman][] | [sjackman.ca][]

[sjackman.github.io/open-science][]

[Fork me on GitHub!][]

[sjackman.github.io/open-science]: https://sjackman.github.io/open-science

Shaun Jackman
------------------------------------------------------------

| [Genome Sciences Centre][], BC Cancer Agency
| Vancouver, Canada
| [\@sjackman][]
| [github.com/sjackman][]
| [sjackman.ca][]

![](images/sjackman.jpg)
