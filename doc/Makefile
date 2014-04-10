html: index.html rdf.html README.html

%.html: %.mdwn pandoc.html
	pandoc --mathjax --template=pandoc.html --css=pandoc.css --toc -o $@ $<

index.html: index.mdwn pandoc.html
	pandoc --mathjax --template=pandoc.html --css=pandoc.css -V pagetitle= -o $@ $<

pandoc.html: contents.mdwn pandoc.html.in
	pandoc --template=pandoc.html.in -o $@ $<

clean:
	$(RM) *.html

.INTERMEDIATE: pandoc.html
