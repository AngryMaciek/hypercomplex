#
# Makefile for installation of the library.
#
# Maciej Bak
# Swiss Institute of Bioinformatics
# 22.10.2019
#

SRC = Hypercomplex.hpp

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    INCLUDE_PREFIX = /usr/include
endif
ifeq ($(UNAME_S),Darwin)
    INCLUDE_PREFIX = /usr/local/include
endif

.PHONY: help install uninstall

# =========================================================
# Print available commands
# =========================================================

help:
	@echo "help - display this message"
	@echo "install - install the library (requires admin rights)"
	@echo "uninstall - uninstall the library (requires admin rights)"

# =========================================================
# Install
# =========================================================

install: $(INCLUDE_PREFIX)/Hypercomplex/$(SRC)

# Create a separate directory for the header-only library
$(INCLUDE_PREFIX)/Hypercomplex:
	mkdir -p $@;

# Copy the header-only library file into the right directory
$(INCLUDE_PREFIX)/Hypercomplex/$(SRC): $(SRC) $(INCLUDE_PREFIX)/Hypercomplex
	cp $< $@

# =========================================================
# Uninstall
# =========================================================

# Remove the library directory
uninstall:
	rm -rf $(INCLUDE_PREFIX)/Hypercomplex
