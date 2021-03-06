# This is a makefile for Peano's StorageSchemes component
# generated by the Peano Development Tools (PDT) 

# Include generated list of compilation targets
-include files.mk

# Minimum java compiler major version requirement for the ExaHyPE toolkit.
# This is Java 8 since September 2017.
EXAHYPE_MINIMUM_JAVA_VERSION := 8

OBJECTS=$(SOURCES:.java=.class)


default: 
	@echo ExaHyPE toolkit makefile
	@echo ------------------------
	@echo 
	@echo Supported targets then are
	@echo   clean
	@echo   createParser
	@echo   compile
	@echo   dist
	@echo   all

# Checks for availability of executables at runtime, without any make target. Source:  https://stackoverflow.com/a/25668869
JAVA_TOOLS = javac javap java jar
JAVA_CHECK_AVAIL := $(foreach exec,$(JAVA_TOOLS),\
  $(if $(shell which $(exec)),[$(exec): yes],$(error "Missing $(exec) in PATH. I need the full JRE and JDK in order to compile and execute the java toolkit. The Path variable is ${PATH}.")))

# This makefile target is called whenever a build fails.
java_version_proof:
	@echo "The build failed and the reason is probably that your Java version is too old. Note that the ExaHyPE Toolkit requires at least Java ${EXAHYPE_MINIMUM_JAVA_VERSION}:"
	@./javac-version-check.sh
	@echo "This was determined on your Path variable ${PATH}"

clean:
	rm -f $(OBJECTS)
	rm -f files.mk
	rm -f ../dist/ExaHyPE.jar
	rm -rf eu/exahype/node

#	delete unpacked dependencies.
	rm -rf org/jtwig
	rm -rf com/github/


createParser:
	java -jar ../lib/sablecc.jar -d . ../exahype.grammar 
	rm -f files.mk


compile: files.mk $(OBJECTS) 

unpack_deps:
#	Ugly trick to get the depencies in the same directory.
#	They only contain .class files, so nothing to compile, needs only
# 	to be present during compiling eu/exahype/.
#
#	jar xf ../lib/jtwig-core.jar org/jtwig
#	jar xf ../lib/mustache.jar   com/github/mustachejava
#	jar xf ../lib/handlebars.jar com/github/jknack

jar:
	jar cmf ./Manifest.txt ExaHyPE.jar .

dist:
	jar cmf ./Manifest.txt ExaHyPE.jar .
	mv ExaHyPE.jar ../dist
	tar -czf ../dist/ExaHyPE-toolkit.tar.gz --exclude=*.class Manifest.txt Makefile eu 


all:
#	try serial making by invoking individual make instead of dependency list:
	make clean
	make createParser || make java_version_proof
	make unpack_deps
	make compile || make java_version_proof
	make dist


files.mk:	
	@touch files.mk
	@echo "SOURCES=\\" > files.mk
	find . -name '*.java' | tr '\n' ' ' >> files.mk



$(OBJECTS): %.class : %.java
	javac -sourcepath . $<  


