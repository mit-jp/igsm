#!/bin/sh
# Assemble examples of ctrparam.sh for a variety of command line permutations to
# configure.

# make sure we're in the top build directory
test -f configure || exit 1
rm -rf ctrparam-test && mkdir ctrparam-test && cd ctrparam-test

# one set of options at a time
cat ../$(dirname $0)/ctrparam.txt | while read opts; do
  # munge the option set for use in a filename
  opt2=`echo $opts | tr -d " " | tr \= -`
  echo -n "Configuring with $opts... "
  # save the configuration output to a .log file
  ../configure $opts >config.$opt2.log
  # copy the generated .h file
  grep "^\#define" src/inc/ctrparam.h >"ctrparam.$opt2.h"
  echo "done"
done
