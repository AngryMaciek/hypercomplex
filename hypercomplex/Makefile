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
	@echo "install - install the library"
	@echo "uninstall - uninstall the library"

# =========================================================
# Install
# =========================================================

install: ~/lib/include/$(SRC)

# Create top-level directory for the libraries
~/lib:
	mkdir -p $@;

# Create top-level directory for the header-only library
~/lib/include:
	mkdir -p $@;

# Copy the header-only library file into the right directory
~/lib/include/$(SRC): $(SRC) ~/lib/include
	cp $< $@

# =========================================================
# Uninstall
# =========================================================

# Remove the top-level directories
uninstall:
	rm -f ~/lib/include/$(HDR);
	if [ $$(ls -a1 ~/lib/include | wc -l)==2 ]; then \
		rm -rf ~/lib/include; \
	fi;
	if [ $$(ls -a1 ~/lib | wc -l)==2 ]; then \
		rm -rf ~/lib; \
	fi;
