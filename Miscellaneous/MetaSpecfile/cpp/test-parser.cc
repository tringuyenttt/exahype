#include <stdio.h>
#include <iostream>

#include "mexa.h"

#include "exahype/Parser.h"

#include <utility> // pair
#include <regex> 

int main(int argc, char** argv) {
	std::string progname = argv[0];
	if (argc < 2) {
		printf("Usage: %s <someSpecfile.exahype>\n", progname.c_str());
		return -1;
	}

	// cmdlineargs contains all argv expect the progname.
	std::vector<std::string> cmdlineargs(argv + 1, argv + argc);
	std::string firstarg = cmdlineargs[0];

	printf("Reading specfile %s...\n", firstarg.c_str());
	
	// run the program with EXAHYPE_VERBOSE_PARSER="TRUE"
	// to see how the parser parses the tokenStream.
	exahype::Parser parser;
	parser.readFile(firstarg);

	if (!parser.isValid()) {
		printf("Invalid specfile\n");
		return -2;
	}
	
	exahype::Parser::ParserView view(parser.getParserView(0));
	
	mexa::mexafile mf = mexa::fromOrderedMap(view.getAllAsOrderedMap(), "specfile");
	
	//std::cout << "mf=" << mf.query_root("boundaries").toString();
	std::cout << "mf=" << mf.toString();
}
