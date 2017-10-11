#!/bin/bash
#
# Outsourced part of the exa.sh: Locating and Determining paths of applications.
#   This used to work based on registered paths where applications are supposed
#   to be stored in directly. We nowadays also support subdirectories.
#

SCRIPT="$(readlink -f $0)" # absolute path to update-installation.sh
ME="$(basename "$SCRIPT")" # just my name (update-installation.sh)
SCRIPTDIR="$(dirname $SCRIPT)"
EXA="$SCRIPTDIR/exa.sh"
GITROOT="$(cd $SCRIPTDIR && git rev-parse --show-toplevel)" # absolute path to ExaHyPE repository working copy

# Paths where applications are typically installed.
# Add your paths white space seperated.
REGISTERED_PATHS="ApplicationExamples Applications AstroApplications Demonstrators"

verbose() { info $@; $@; } # only used for debugging
info () { echo -e $ME: $@; } # print error/info message with name of script
fail () { info $@; exit -1; } # exit after errormessage with name of script
abort () { echo -e $@; exit -1; } # fail without name of script
finish () { echo $@; exit 0; } # finish with message happily
subreq() { $SCRIPT $@; } # subrequest: Query another command for output
cdroot() { cd "$GITROOT"; } # the crucial change to the repository root directory

CMD="$1" # the actual command
PAR="$2" # some parameter (for passing to bash functions)

# note that we have "set -e" here in contrast to exa.sh. Therefore
# make sure that every expression return well-understood.

# getappname: set $APPNAME or die
getappname() {
	# Try to read the application name from the argument.
	APPNAME="$PAR";
	[[ -z "$APPNAME" ]] && abort "Usage: $0 $CMD <AppName>" || true;
	# Sanatize it: Remove the .exahype extension, if exists
	APPNAME="${APPNAME%.exahype}"
}

# test for files/directories and print + exit or return.
testfile() { ls -f "$@" 2>/dev/null && exit 0 || true; }
testdir()  { ls -d "$@" 2>/dev/null && exit 0 || true; }

originalpwd="$PWD"; cdroot
set -e

case $CMD in
	"tree"|"all") # List all ExaHyPE applications in a nice way.
		info "Listing available Applications:"
		if which tree 2>&1; then
			tree -d -L 1 $REGISTERED_PATHS
		else
			# poor man's tree. Does not work 100%
			for path in $REGISTERED_PATHS; do
				echo "*** $path ***"
				[[ -e $path ]] && find $path -type d -exec basename {} \; || info "Does not exist. Please remove from $ME"
				
			done
		fi

		info "Showing also all specfiles below ExaHyPE's root:"
		find . -type f -iname '*.exahype'
		;;
	"list") # Lists all ExaHyPE applications machine readable. Use "find" for full path.
		for path in $REGISTERED_PATHS; do
			[[ -e $path ]] && for appdir in $(ls -d $path/*/); do
				echo "$appdir"
			done
		done
		;;
	"list-abs") # List with absolute paths. Yeah, this should be an argument --abs instead.
		# Option 1 (gives paths like /....foo/bar/)
		# prepend="$GITROOT/"
		# subreq list | awk "{ print \"$prepend\" \$0 }"
		
		# Option 2: (gives canonical paths like /...foo/bar)
		for path in $REGISTERED_PATHS; do
			[[ -e $path ]] && for appdir in $(ls -d $path/*/); do
				realpath "$appdir"
			done
		done
		;;
	"app"|"appdir") # Gives the full path from ExaHyPE root to an application
		getappname
		specfile=$(subreq specfile $APPNAME)

		# This is a naive parser. If it fails too often, one should instead rely on the proper
		# C++ or Java specfile parser. It is probably even useful to write an own standalone
		# C++ parser which converts the specfile into JSON which can be queried with jq.
		
		cat $specfile | awk 'BEGIN{r=1} /output-directory/{ r=0; print $4; } END{ exit r}' \
			|| fail "Could not extract application directory from specfile '$specfile' which I identified for application '$APPNAME'"
		;;
	"spec"|"specfile") # Gives the full path from ExaHyPE root to an application specfile
		getappname
		for path in $REGISTERED_PATHS; do
			# allow several variants:
			testfile $path/$APPNAME.exahype
			testfile $path/$APPNAME/$APPNAME.exahype
			testfile $path/*/$APPNAME.exahype
			testfile $path/*/*/$APPNAME.exahype
		done
		fail "Application '$APPNAME' not found";
		;;
	"fuzzy-specfile") # List all specfiles which match the asked one. Don't rely on this in scripts.
		getappname
		# First, ask for a specific specfile. If one matches, this is good.
		
		# Second, look for candidates. Echo them all. Do not exit after individual ones:
		for path in $REGISTERED_PATHS; do
			# collect all candidates
			echo TODO
		done
		;;
	"fuzzy-app") # The same as fuzzy-specfile, just printing app directories
		echo TODO
		;;
	"parent") # Gives directory where app lives inside.
		# This command used to be "find-appdir" before!! It was renamed to "parent" because
		# the actual parent dir of an application no more guarantees that the specfile is in the same directory.
		getappname
		for path in $REGISTERED_PATHS; do
			ls -d "$path/$APPNAME" &>/dev/null && finish "$path/" || continue
		done
		fail "Could not find Application '$APPNAME' somewhere"
		;;
	"binary") # Gives path to the executable, even if not present, by inspecting the specfile.
		getappname;
		SPECFILE="$(subreq specfile "$APPNAME")" || abort "Specfile Failure: $SPECFILE"
		PROJECTNAME=$(grep '^exahype-project' ${SPECFILE} | awk '{ print $2; }')
		APPPATH="$(subreq app "$APPNAME")"
		echo $APPPATH/ExaHyPE-$PROJECTNAME
		;;
	"reverse") # Determine name of application based on argument path
		TESTDIR="$PAR"; [[ -z "$TESTDIR" ]] && abort "Usage: $0 reverse <SomePath>"
		# This is a stupid search. It will work only with absolute canonical paths,
		# ie. cannot detect symlinks. You should pass the argument throught `readlink`
		# before, cf. the reverse-pwd subcommand.
		
		# @TODO: This is still open and does not work properly. ALERT
		
		subreq list-abs | while read appline; do
			# fun thing: echo -e "1\n2\n3" | while read foo; do echo $foo; exit; done;
			# this does not exit the script but only the while loop. This is weird.
			grep -q -i "$appline" <<< "$TESTDIR" && { echo "$(basename "$appline")"; exit -1; }
		done
		fail "Could not find $TESTDIR in the applications, cf. output of 'exa find list-abs' for reference list"
		;;
	"reverse-pwd") # Determine the name of an application based on the current working directory
		# bring to canonical form:
		pwd="$(realpath "$originalpwd")"
		info "Looking for application residing in $pwd"
		exec $SCRIPT reverse "$pwd" # no subreq for exec
		;;
	""|"help"|"--help")
		info "Use this script to locate applications. The available commands are:"
		info "    list, appdir, specfile, parent, binary"
		info "Try '$ME list' to see all applications available."
		info "Try '$ME appdir EulerFlow' to learn the path to the EulerFlow application."
		info "Try '$ME spec EulerFlow' to learn the path to the specfile of EulerFlow"
		info "etc."
		# todo: Write better help.
		;;
	*)
		fail "Could not understand command '$CMD'"
		;;
esac
