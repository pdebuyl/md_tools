html: index.html rdf.html README.html

MATHJAX_URL=/usr/share/javascript/mathjax/unpacked/MathJax.js

%.html: %.mdwn pandoc.html
	pandoc --mathjax --template=pandoc.html --css=pandoc.css --toc -o $@ $<

index.html: index.mdwn pandoc.html
	pandoc --mathjax=$(MATHJAX_URL) --template=pandoc.html --css=pandoc.css -V pagetitle= -o $@ $<

pandoc.html: contents.mdwn pandoc.html.in
	pandoc --template=pandoc.html.in -o $@ $<

clean:
	$(RM) *.html

.INTERMEDIATE: pandoc.html
