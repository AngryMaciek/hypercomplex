#!/bin/bash

# create a non-root user "user" passed from the command line (or fallback);
# exec the rest (CMD) as "user"
# ~AngryMaciek

ID=${HOSTUID:-9001}
useradd --shell /bin/bash -u $ID -o -c "" -m user
export HOME=/home/user
exec /usr/sbin/gosu user "$@"
