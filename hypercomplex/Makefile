#
# Makefile for installation & development of the library.
#
# Maciej Bak
# Department of Mathematics, City University of London
# 22.10.2019
#

SRC1 = Hypercomplex.hpp
SRC2 = Polynomial.hpp

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    INCLUDE_PREFIX = /usr/include
endif
ifeq ($(UNAME_S),Darwin)
    INCLUDE_PREFIX = /usr/local/include
endif

.PHONY: help install uninstall test lint docs

# =========================================================
# Print available commands
# =========================================================

help:
	@echo "help - display this message"
	@echo "install - install the library (requires admin rights)"
	@echo "uninstall - uninstall the library (requires admin rights)"
	@echo "test - execute testing framework (requires mpfr library)"
	@echo "lint - run static code analysis (requires cpplint)"
	@echo "docs - generate project's documentation (requires doxygen)"

# =========================================================
# Install
# =========================================================

install: $(INCLUDE_PREFIX)/Hypercomplex/$(SRC1) $(INCLUDE_PREFIX)/Hypercomplex/$(SRC2)

# Create a separate directory for the header-only library
$(INCLUDE_PREFIX)/Hypercomplex:
	mkdir -p $@;

# Copy the header-only library file into the right directory
$(INCLUDE_PREFIX)/Hypercomplex/$(SRC1): $(SRC1) $(INCLUDE_PREFIX)/Hypercomplex
	cp $< $@

# Copy the helper library file into the right directory
$(INCLUDE_PREFIX)/Hypercomplex/$(SRC2): $(SRC2) $(INCLUDE_PREFIX)/Hypercomplex
	cp $< $@

# =========================================================
# Uninstall
# =========================================================

# Remove the library directory
uninstall:
	rm -rf $(INCLUDE_PREFIX)/Hypercomplex

# =========================================================
# Test
# =========================================================

# Prepare, compile, test, cleanup
test:
	mkdir ../.test/unit/hypercomplex; \
	cp Hypercomplex.hpp ../.test/unit/hypercomplex/Hypercomplex.hpp; \
	cp Polynomial.hpp ../.test/unit/hypercomplex/Polynomial.hpp; \
	cd ../.test/unit; \
	g++ -O2 -Wall --std=c++17 -o test test.cpp -lmpfr -lgmp; \
	./test -d yes -w NoAssertions --use-colour yes --benchmark-samples 100 --benchmark-resamples 100000; \
	rm -rf hypercomplex test

# =========================================================
# Lint
# =========================================================

# Run static code analysis
lint:
	cpplint Hypercomplex.hpp
	cpplint Polynomial.hpp

# =========================================================
# Docs
# =========================================================

# Generate doxygen documentation 
docs:
	cd ..; \
	doxygen Doxyfile
