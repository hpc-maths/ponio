# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS  =
SPHINXBUILD = sphinx-build
SOURCEDIR   = source
BUILDDIR    = build
SOURCERST   = $(shell find source -name *.rst) 

include Doxyfile

all: html

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile $(XML_OUTPUT)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
$(XML_OUTPUT): $(wildcard $(INPUT)/$(PROJECT_NAME)/*.hpp)
	doxygen

clean:
	rm -rf $(XML_OUTPUT) $(BUILDDIR)

html: $(XML_OUTPUT) $(SOURCERST)
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

