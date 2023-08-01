# Based on:
# https://www.gitpod.io/docs/configure/workspaces/workspace-image#custom-base-image
FROM ubuntu:22.04

# Install:
# - git (and git-lfs), for git operations (to e.g. push your work).
#   Also required for setting up your configured dotfiles in the workspace.
# - gosu, while not required, is recommended to be installed, since the
#   workspace user (`gitpod`) is non-root and won't be able to install
#   and use `sudo` to install any other tools in a live workspace.
RUN apt-get update && apt-get install -yq --no-install-recommends \
    git git-lfs gosu g++ make valgrind libmpfr-dev python3 python3-pip doxygen \
    && pip install cpplint pre-commit \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/*

# Create the gitpod user. UID must be 33333.
RUN useradd -l -u 33333 -md /home/gitpod -s /bin/bash -p gitpod gitpod

USER gitpod
