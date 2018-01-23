#!/bin/bash
#
# Use this script as a loose coupling of the ExaHyPE-Engine git
# repository to the Peano svn repository.
#

# cd to the Peano basedir. Thus script can be executed from anywhere.
cd $(dirname "$0")

SVN_URL="svn://svn.code.sf.net/p/peano/code/trunk/src"

# if you experience the svn protocol blocked or not working, try this.
# Source: https://sourceforge.net/p/peano/code/HEAD/tree/
SVN_ALT_URL="https://svn.code.sf.net/p/peano/code/trunk/src"

if test -e .svn
then
	echo "Updating peano"
	exec svn update
else
	echo "Checking out peano from $SVN_URL"
	if ! svn checkout $SVN_URL .
	then
		echo "URL not working, trying the alternative one $SVN_ALT_URL"
		exec svn checkout $SVN_ALT_URL .
	fi
fi
