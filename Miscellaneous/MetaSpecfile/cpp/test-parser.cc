#include <stdio.h>

// back to access private members for the time being.
//#define private public
//#define protected public

#include "exahype/Parser.h"

#include <utility> // pair

typedef std::pair<std::string, std::string> kv;
typedef std::vector<kv> ordered_map; // poor man's.

/*
ordered_map get_all_tokens(exhype::Parser& _parser, const std::string& key) const {
  std::string token;
  std::regex  COLON_SEPARATED(R"((.+):(.+))");
  std::smatch match;

  token = _parser.getTokenAfter("solver", _solverNumberInSpecificationFile + 1, "constants", 1, 0);
  std::regex_search(token, match, COLON_SEPARATED);

  int i = 1;
  while (match.size() > 1) {
    if (match.str(1).compare(key)==0) {
      logDebug("getValue()", "solver " << _solverNumberInSpecificationFile + 1 << ": found constant '" << key << "' with value '" << match.str(2) << "'.");
      return match.str(2);
    }
    token = _parser.getTokenAfter("solver", _solverNumberInSpecificationFile + 1, "constants", 1, i++);
    std::regex_search(token, match, COLON_SEPARATED);
  }
  logDebug("hasKey()", "solver " << _solverNumberInSpecificationFile + 1 << ": cannot find constant '" << key << "'.");
  return "";
}
*/

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

}
