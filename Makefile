###############################################################################
#
#   Makefile for installation & development of the library.
#
#   AUTHOR: Maciek Bak
#   AFFILIATION: Department of Mathematics, City University of London
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 22.10.2019
#
###############################################################################

SRC1 = Hypercomplex.hpp
SRC2 = Polynomial.hpp

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    INCLUDE_PREFIX = /usr/include
endif
ifeq ($(UNAME_S),Darwin)
    INCLUDE_PREFIX = /usr/local/include
endif

.PHONY: help install uninstall test lint docs build clean image container

# =============================================================================
# Print available commands
# =============================================================================

help:
	@echo "help - display this message"
	@echo "install - install the library (requires admin rights)"
	@echo "uninstall - uninstall the library (requires admin rights)"
	@echo "test - execute testing framework (requires mpfr library)"
	@echo "lint - run static code analysis (requires cpplint)"
	@echo "docs - generate project's documentation (requires doxygen)"
	@echo "build - build conda package (requires conda-build)"
	@echo "clean - remove all dev artifacts"
	@echo "image - build Docker image (requires Docker deamon running)"
	@echo "container - start interactive Docker container (requires image)"

# =============================================================================
# Install
# =============================================================================

install: $(INCLUDE_PREFIX)/Hypercomplex/$(SRC1) $(INCLUDE_PREFIX)/Hypercomplex/$(SRC2)

# Create a separate directory for the header-only library
$(INCLUDE_PREFIX)/Hypercomplex:
	@mkdir -p $@

# Copy the header-only library file into the right directory
$(INCLUDE_PREFIX)/Hypercomplex/$(SRC1): hypercomplex/$(SRC1) $(INCLUDE_PREFIX)/Hypercomplex
	@cp $< $@

# Copy the helper library file into the right directory
$(INCLUDE_PREFIX)/Hypercomplex/$(SRC2): hypercomplex/$(SRC2) $(INCLUDE_PREFIX)/Hypercomplex
	@cp $< $@

# =============================================================================
# Uninstall
# =============================================================================

# Remove whole library directory
uninstall:
	@rm -rf $(INCLUDE_PREFIX)/Hypercomplex

# =============================================================================
# Test
# =============================================================================

# Prepare, compile, test, cleanup
test:
	@mkdir .test/unit/hypercomplex
	@cp hypercomplex/Hypercomplex.hpp .test/unit/hypercomplex/Hypercomplex.hpp
	@cp hypercomplex/Polynomial.hpp .test/unit/hypercomplex/Polynomial.hpp
	@g++ -O0 -Wall --std=c++17 -o test .test/unit/test.cpp -lmpfr -lgmp
	@./test [unit] -d yes -w NoAssertions --use-colour yes --benchmark-samples 100 --benchmark-resamples 100000
	@rm -rf .test/unit/hypercomplex test

# =============================================================================
# Lint
# =============================================================================

# Run static code analysis
lint:
	@cpplint hypercomplex/Hypercomplex.hpp
	@cpplint hypercomplex/Polynomial.hpp

# =============================================================================
# Docs
# =============================================================================

# Generate doxygen documentation
docs:
	@doxygen Doxyfile

# =============================================================================
# Build
# =============================================================================

# Build conda package
build:
	@conda build . -c conda-forge

# =============================================================================
# Clean
# =============================================================================

# Remove all artifacts
clean:
	@find . -type f -name '*.DS_Store' -delete
	@rm -rf docs/html docs/latex
	@rm -f test

# =============================================================================
# Image
# =============================================================================

# Build Docker image
image:
	@# Build docker image
	@if [ $$(docker ps &> /dev/null; echo $$?) -eq 0 ]; then\
		docker rm hypercomplex &> /dev/null;\
		sleep 1;\
		docker rmi hypercomplex:latest &> /dev/null;\
		docker build -t hypercomplex:latest -f Dockerfile .;\
	else \
		echo "Docker deamon is not running OR the current user requires higher privileges.";\
	fi

# =============================================================================
# Container
# =============================================================================

# start Docker container
container:
	@# Run a docker container
	@if [ $$((docker images | grep hypercomplex) &> /dev/null; echo $$?) -eq 0 ] \
	; then\
		docker rm hypercomplex &> /dev/null;\
		docker run --name hypercomplex -e HOSTUID=`id -u $$USER` -p 8888:8888 -it -v $$PWD:/hypercomplex hypercomplex:latest;\
	else \
		echo "Docker deamon isn't running OR hypercomplex image isn't available OR the current user requires higher privileges.";\
	fi
