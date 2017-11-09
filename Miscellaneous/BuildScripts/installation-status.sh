#!/bin/bash
#
# Dump status information about changes in the current ExaHyPE installation.
# This can be helpful for debugging. Output is made to stdout. No files are
# created.
#
# This script looks both into it's own surrounding area as well as the current
# working directory.
#
# Written by SvenK at 2017-11-06.
#

BUILDSCRIPTS="$(dirname "$0")"

have() { which $@ 2>/dev/null >/dev/null; }

if have git; then
	echo "$(git --version) is installed on $(hostname)."

	# first, check for local changes
	if [ -d .git ] || git rev-parse --git-dir > /dev/null 2>&1; then 
		echo "The current working directory $PWD is under GIT version control"
		echo "at HEAD commit $(git show --oneline -s)"
		echo "Local directory changes:"
		git status .
	fi

	# second, check the buildscripts
	cd "$BUILDSCRIPTS"
	if [ -d .git ] || git rev-parse --git-dir > /dev/null 2>&1; then 
		GITROOT="$(git rev-parse --show-toplevel)"
		echo "Found git repository location at $GITROOT"
		echo "Current HEAD commit:   $(git show --oneline -s)"
		echo "On branch:"
		git branch
		echo "Changes made at the global repository:"
		cd $GITROOT
		git status
	else
		echo "Could not find a git repository at the Buildscripts $BUILDSCRIPTS"
	fi
else
	echo "No 'git' on the PATH. Cannot do any statement."
fi

if have svn; then
	echo "$(svn --version | head -n1) is installed on $(hostname)"
	echo "@TODO: Check and dump local Peano changes." # TODO
fi
