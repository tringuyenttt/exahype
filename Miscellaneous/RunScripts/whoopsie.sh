#!/bin/bash
#
## This is a more generic version of Misc/BuildScripts/whoopsie-paster.sh, which is supposed
## to dump all kind of debug-related data. Instead, this script tries to do so in terms 
## of ExaHyPE executables.
##
## The short idea of this script is to dump all kind of information in a text file and upload
## it to a pastebin service such as sprunge.us or termbin.com. If this is not possible, instead
## provide the user the generated file for dumping on his own.
##
## Usage examples:
##
##    ../path/to/whoopsie.sh  ./ExaHyPE-Euler ../EulerFlow.exahype
##    exa whoopsie ./ExaHyPE-Euler ../EulerFlow.exahype
##    exa whoopsie mpirun -np 10 ./ExaHyPE-Euler ../EulerFlow.exahype
##
## Minimalistic demonstrators for showing how success and failure is captured:
##
##    exa whoopsie true
##    exa whoopsie false
##
## Written by SvenK at 2017-11-06.
#

dumplog="dump.log"
# wipe the file
echo -n > $dumplog

# self documenting: if no parameter given
if [ $# -eq 0 ]; then grep -E '^##' "$0"; exit; fi

have() { which $@ 2>/dev/null >/dev/null; } # checks on the path, not bash builtins

runscripts="$(dirname "$0")"
buildscripts="${runscripts}/../BuildScripts"
exa="${buildscripts}/exa.sh"

# Turn off buffering in pipes. Set this to the empty string if your
# system does not have the command "stdbuf"
if have stdbuf; then unbuf="stdbuf -i0 -o0 -e0"; else unbuf=""; fi

# Use like "something | $teelog": Dumps both into log and stdout.
teelog="$unbuf tee -a $dumplog"
hiddenexec() { $@ >> $dumplog 2>&1; } # execute something and put all output into dumplog
log() { hiddenexec echo $@; } # log something
logcmd() { log ">>> $@"; } # log a command
verbose() { logcmd $@; hiddenexec $@; } # verbose a command and execute it
log_file() { log "----[ begin of $@ ]----"; hiddenexec cat $@; log "----[ end of $@ ]----"; } # log a text file
spacing() { log; } # just some lines of spacing


log "This is the ExaHyPE 'whoopsie' runtime problem catcher. I accompany"
log "the user $(whoami) on the computer node $(hostname) at"
log "time $(date) by running ExaHyPE. Running ExaHyPE is not easy."
spacing
log "So the user starts with an environment given by"
spacing
verbose $exa check
spacing
log "You also might to want to store coredumps:"
verbose ulimit -c   # ulimit is bash builtin
spacing
log "Now try to run the command:"
spacing
spacing
# Try to run the commands, capture also the output. If it finishes: fine.
logcmd $@
set -o pipefail  # needed for detecting failure in pipes
if 2>&1 $@ | $teelog; then
	log "Finished successfully."
	exit
fi
# else: Log all possible stuff
echo "Whoopsie: Catching a failed command to $PWD/$dumplog"
log ">>> The command FAILED. " # with return value $? <= not the return value of the program. Needs Pipe inspection.
spacing
log "Post mortem inspection of the command '$@' by looking into each of its words:"
spacing

for part in $@; do
	if [[ -e $part ]]; then
		if [[ -x "$part" ]] && ./$part --help 2>&1 | grep -qi exahype; then
			log "=> $part is an ExaHyPE executable. This is how it was compiled:"
			log
			verbose $part --version
		elif have file; then
			if file $part | grep -qi ascii; then
				log "=> $part is a text file. Here are it's contents:"
				log_file $part
			else
				log "=> $part is a file but I don't know what's inside. This is what I can learn about it:"
				verbose file $part
			fi
		else
			log "=> $part is a file but I miss the Unix tool 'file' to look into it."
		fi
	elif have which; then
		# only check parts which don't go like "-foo" or "--foo" because which interpretes this as parameters
		if ! [[ $part == -* ]] && which $part; then
			log "=> $part is on the PATH. It resolves to:"
			verbose which $part
		else
			log "=> '$part' does not seem to be on the PATH nor it is a file. Probably a plain old parameter."
		fi
	else
		log "=> '$part' is not a file and I have no tools (no 'which') to find out what it is."
	fi
done

spacing
spacing
log "The current environment: "
spacing
verbose env
spacing
log "Inspection of local changes in the installation:"
verbose $buildscripts/installation-status.sh

log
if [[ -e make.log ]]; then
	log "I have found a make.log which probably explains how this was built:"
	log_file make.log
else
	log "No make log found. Are there other log files? Let's check:"
	verbose ls *.log
fi

uploadWith() {
	echo "Uploading the file ${dumplog}, please send the following link to your collaborators:"
	cat $dumplog | $@ || { echo "Failure with uploading! Please just send the file $PWD/${extendedlog} as an email attachment to your collaborators."; }
}

if have curl; then
	uploadWith 'curl -s -F sprunge=<- http://sprunge.us'
elif have nc; then
	uploadWith 'nc termbin.com 9999'
else
	echo "Please send the file $PWD/${dumplog} as an email attachment to your collaborators."
fi

