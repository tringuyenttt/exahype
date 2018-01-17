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
	
	mexa::mexafile mf = mexa::fromSpecfile(view.getAllAsOrderedMap(), "specfile");
	
	std::cout << "mf=" << mf.toString();
	
	// Testing queries
	mexa::mexafile mq = mf.query_root("boundaries");
	std::cout << "mq=" << mq.toString() << "\n";
	std::cout << "contains left? " <<  mq.contains("left") << "\n";
	std::cout << "value from left = <" << mq.get("left").get_string() << ">\n";
	
	// testing vector
	/*
	std::vector<double> pos = mf.vec("shocktube/right/bmag", 3).get_double();
	//std::vector<int> pos = mf.vec("shocktube/right/bmag", 3).as_int();
	for(auto j : pos)
		std::cout << "vector value: " << j << "\n";
	*/
	
	// Testing strings
	/*
	mexa::mexafile mu  = mf.query_root_require("initialdata");
	std::cout << "mu=" << mu.toString();
	std::cout << "named=" << mu("name").toString();
	*/
}
