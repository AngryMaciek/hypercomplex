#
# Makefile for build and
# installation of the C++ library.
#
# Maciej Bak
# Swiss Institute of Bioinformatics
# 22.10.2019
#

CC = g++
CFLAGS = -O2
SRC = Hypercomplex.cpp
HDR = $(SRC:.cpp=.h)
OBJ = $(SRC:.cpp=.o)
LIB = $(OBJ:.o=.a)

.PHONY: all install clean uninstall

# =========================================================
# Build
# =========================================================

all: $(LIB)

# Create an object file
$(OBJ):
	$(CC) $(CFLAGS) -c -o $@ $(SRC)

# Create a static library file
$(LIB): $(OBJ)
	ar rcs $@ $<

# =========================================================
# Install
# =========================================================

install: ~/lib/lib/$(LIB) ~/lib/include/$(HDR)

# Create top-level directory for the libraries
~/lib:
	mkdir -p $@;

# Create top-level directory for the static library files
~/lib/lib:
	mkdir -p $@;

# Create top-level directory for the headers
~/lib/include:
	mkdir -p $@;

# Copy the library file into the right directory
~/lib/lib/$(LIB): $(LIB) ~/lib/lib
	cp $< $@

# Copy the header file into the right directory
~/lib/include/$(HDR): $(HDR) ~/lib/include
	cp $< $@

# =========================================================
# Clean
# =========================================================

# Remove the build files
clean:
	rm -f $(OBJ);
	rm -f $(LIB);

# =========================================================
# Uninstall
# =========================================================

# Remove the top-level directories
uninstall:
	rm -f ~/lib/lib/$(LIB);
	rm -f ~/lib/include/$(HDR);
	if [ $$(ls -a1 ~/lib/lib | wc -l)==2 ]; then \
		rm -rf ~/lib/lib; \
	fi;
	if [ $$(ls -a1 ~/lib/include | wc -l)==2 ]; then \
		rm -rf ~/lib/include; \
	fi;
	if [ $$(ls -a1 ~/lib | wc -l)==2 ]; then \
		rm -rf ~/lib; \
	fi;
