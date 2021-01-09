#
# Makefile for installation of the library.
#
# Maciej Bak
# Swiss Institute of Bioinformatics
# 22.10.2019
#

SRC = Hypercomplex.h

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

install: /usr/include/Hypercomplex/$(SRC)

# Create a separate directory for the header-only library
/usr/include/Hypercomplex:
	mkdir -p $@;

# Copy the header-only library file into the right directory
/usr/include/Hypercomplex/$(SRC): $(SRC) /usr/include/Hypercomplex
	cp $< $@

# =========================================================
# Uninstall
# =========================================================

# Remove the library directory
uninstall:
	rm -rf /usr/include/Hypercomplex
