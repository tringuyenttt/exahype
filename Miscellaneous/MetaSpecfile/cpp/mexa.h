#ifndef MEXA_META_EXAHYPE_PARFILE_CPP_IMPL_HEADER
#define MEXA_META_EXAHYPE_PARFILE_CPP_IMPL_HEADER

// C++ STL:
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <map>
#include <utility>
#include <vector>
#include <cstring>
#include <string>
#include <cstdarg>

/**
 * Mexa is a C++11 standalone implementation for the simple-mexa format. This format
 * consists of line-by-line text files in the format
 *
 *   foo/baz = "bar"
 *   foo/bar = 0.9
 *   foo/foo = Yes
 *
 * i.e. it knows hierarchy, assignments and datatypes. In contrast to the fully
 * specified mexa format (which understands inclusion, extension, symbolic variables,
 * abstract trees/lists), this format is very simple.
 * For details about the file format, read the mexa README, for instance at
 * https://bitbucket.org/svek/mexa
 *
 * Written by SvenK in Nov 2017 for ExaHyPE.
 **/
namespace mexa {
	typedef std::vector<std::string> stringvec;
	struct mexafile; // forward definition

	/**
	 * A symbol represents an hierarchic path. Such a path can point either
	 * to a node in the file or to a leaf where actual information are stored.
	 * This class allows several arithmetic operations with these paths.
	 **/
	struct symbol {
		typedef stringvec path_t;
		path_t path;
		static path_t parseSymbol(std::string name);
		
		symbol() {} ///< the root symbol
		symbol(const std::string& name) : path(parseSymbol(name)) {}
		symbol(const char* const name) : symbol(std::string(name)) {} // for convenience
		symbol(const path_t& path) : path(path) {}
		
		/// Canonical representation of symbol
		std::string toString() const;
		
		/// Lists all parenting symbols, for instance for /a/b/c it is the list [/, /a, /a/b].
		std::vector<symbol> ancestors() const;
		
		/// Returns a new symbol which is prefixed by this
		symbol prefix_add(const symbol& other) const;
		
		/// Remove common prefix
		symbol prefix_remove(const symbol& other) const;
		
		/// For the usage as key in std::map
		bool operator <(const symbol& rhs) const;
		
		/// For comparison of equality (contains, ==)
		bool operator==(const symbol& rhs) const;
	};
	
	
	/// The sourcemap represents a line in a file for understandable error messages.
	struct sourcemap {
		std::string filename;
		int linenumber;
		// Quickly create an "unknown" file
		sourcemap() : filename("unknown"), linenumber(0) {}
		sourcemap(std::string filename, int linenumber) : filename(filename), linenumber(linenumber) {}
		std::string toString() const;
	};
	
	/**
	 * The value class holds either an int, double, bool, string or an undefined
	 * "null"-type value. I did not use
	 *   - Unions because they only store PODs, so the std::string would have
	 *     to be saved out of the union.
	 *   - std::variant<int,double,bool,string> because it is C++17 and I want
	 *     this code to be C+11 compatible.
	 *   - Any compromise with unions because I don't really care about memory
	 *     or speed in this mexa implementation.
	 * 
	 * Furthermore to the union job, this class also allows casting of values.
	 * This may look unneccessary in scalar context, i.e you can always write
	 * 
	 *   double foo = value(3.14).get_int()
	 * 
	 * but it can turn out handy in the vectorial context when a vector consists
	 * of values such as [ 0, 0, 4.5 ], i.e. [int, int, double] but you want
	 * a double vector. See the vector_value class below for an example
	 * usage.
	 **/
	class value {
	public:
		enum class Type { INT, DOUBLE, BOOL, STRING, UNDEF };
	private:
		// if you want, join these to a union or so.
		int i;
		double d;
		bool b;
		std::string s;
		Type a;
	public:
		value(const value&) = default;
		value& operator=(const value&) = default;
		
		value() : a(Type::UNDEF) {}
		value(int i) : i(i), a(Type::INT) {}
		value(double d) : d(d), a(Type::DOUBLE) {}
		value(bool b) : b(b), a(Type::BOOL) {}
		value(std::string s) : s(s), a(Type::STRING) {}
	
		/// Returns whether the instance holds a given type.
		bool isActive(Type type) const { return a == type; }
		/// Checks if this instance holds a given type. If not, raises an exception.
		void assertActive(Type type) const;
		
		// The get_* methods only work if the type is *exactly*
		// what you asked for. There is no casting taking place, i.e.
		// you cannot get_double() for something which is an integer.
		
		int get_int() const { assertActive(Type::INT); return i; }
		double get_double() const { assertActive(Type::DOUBLE); return d; }
		bool get_bool() const { assertActive(Type::BOOL); return b; }
		std::string get_string() const { assertActive(Type::STRING); return s; }
		
		// Instead, the as_* methods allow casting. This works between
		// the {int,double,bool} types. In contrast, the string type
		// is always a string representation of the value. However,
		// a string will never be casted to any of int,bool,double.
		// This is the job of the parser when creating a mexafile.
		
		int as_int() const;
		double as_double() const;
		bool as_bool() const;
		std::string as_string() const;
		
		/// To support this idiom, the canCast function explains what can
		/// be casted and what not. To be used as a replacement for isActive().
		bool canCastTo(Type type) const;
		
		/// Internally, this method does the job. The type information
		/// is only needed for debugging purposes.
		template<typename T> T cast_as_numerical(Type type) const;
		
		
		/// A string representation of the content.
		/// Will read like "bool(True)" or "string(Bla)"
		/// Do not mix this up with get_string() which gives you the
		/// actual content, i.e. "Bla" in the given example (without quotes).
		std::string toString() const;
		
		/// Jus as a helper: A type to string
		static std::string type2str(Type type);
	};
	
	/// a vector value class, forward definition
	class vector_value;
	
	/**
	 * The assignment is the simple and only operation of the simple mexa file
	 * format. It has a LHS (the symbol key) and RHS (a primitive value). It
	 * also holds a reference to the source to allow printing useful error messages
	 * for the user.
	 * This is a Plain old datatype (POD).
	 **/
	struct assignment {
		symbol key;
		value val;
		sourcemap src;
		
		/// A constructor like function to preserve the POD status.
		static assignment make(symbol key, value val, sourcemap src);
		std::string toString() const;
	};
	
	/**
	 * The mexafile represents a parameter file in the simple-mexa file format, i.e. it
	 * is only a map of symbols onto lines/strings (with source information for speaking
	 * error messages).
	 *
	 * This is a plain old datatype, it can be constructed empty. There are functions
	 * which return this structure, parsed from a file.
	 *
	 * The class supports various querying mechanisms inspired by XML-Xpath as well as
	 * actual getters for several C++ types. In case of errors, we rely on exceptions.
	 **/
	struct mexafile {
		/**
		 * The list of key-value pairs is our implementation of an Ordered Map
		 * which replaces the more obvious std::map<symbol, sourced>.
		 **/
		std::vector<assignment> assignments;
		
		/**
		 * Add an item to the ordered assignment list.
		 **/
		void add(const symbol key, const value val, const sourcemap src);
		
		/**
		 * Obtain an item from the ordered assignment list.
		 * We return a copy in order to let this be a const function.
		 **/
		value get(const symbol key) const;
		
		/// Allow also access with call semantics
		value operator()(const symbol key) const { return get(key); }
		/// Allow also access with vector semantics
		value operator[](const symbol key) const { return get(key); }
		
		// High level access:
		vector_value vec(const symbol node, size_t required_length=-1) const;
		vector_value vec(size_t required_length=-1) const;
				
		/**
		 * Returns a string representation which actually is a valid simple-mexa
		 * file itself (if the input is correct).
		 **/
		std::string toString() const;

		/**
		 * Returns true if this file holds the questioned symbol.
		 * @arg doRaise if true, raise an exception if symbol not given
		 **/
		bool contains(symbol leaf, bool doRaise=false) const;
		
		/**
		 * Returns true if the assignment list is empty.
		 **/
		bool isEmpty() const { return assignments.empty(); }

		///////////////////////////////
		// querying
		///////////////////////////////
		
		/**
		 * Filter the file for a given symbol. Returns a new mexafile.
		 **/
		mexafile query(const symbol root) const;
		
		/**
		* Query and give relative paths to root. Ie. if assignments="a/x=1,a/y=2,b/z=3",
		* then query(a) = "a/x=1,a/y=2" while query_root(a)="x=1,y=2".
		**/
		mexafile query_root(const symbol root) const;
		
		/**
		 * As query_root, but aborts if nothing found below the symbol.
		 **/
		mexafile query_root_require(const symbol root) const;

		/**
		 * Returns a copy of this mexafile where the symbol `root' is added to
		 * every entry in the file.
		 **/
		mexafile prefix_add(const symbol root) const;

		/**
		 * Returns a copy of this mexafile where the symbol `root' is subtracted
		 * from every entry in the file.
		 **/
		mexafile prefix_remove(const symbol root) const;
	};
	
	/**
	 * This class mimics the behaviour of the value class. However, it will return
	 * STL vector instances instead. Instances of this class are created by the
	 * mf::vec() method.
	 **/	
	class vector_value {
		mexafile mf, mq;
		symbol node;
		size_t required_length;
		
		/**
		 * Calls the get_*() method from the value class on every value
		 * in the mf.query(node).assignments list. Checks if all types
		 * match, otherwise raises an exception.
		 **/
		template<typename T>
		std::vector<T> get( T (value::*getter)() const, value::Type type, bool doCast) const;
	public:
		vector_value(mexafile mf, symbol node, size_t required_length);
		
		// if you want to be a harsh vector and have all entries to be the same type
		std::vector<int> get_int() const;
		std::vector<double> get_double() const;
		std::vector<bool> get_bool() const;
		std::vector<std::string> get_string() const;
		
		// if instead you allow casting, i.e. in a vector with [0,0,4.5] you
		// could read this as ints to [0,0,4] or as doubles to [0.0,0.0,4.5].
		std::vector<int> as_int() const;
		std::vector<double> as_double() const;
		std::vector<bool> as_bool() const;
		std::vector<std::string> as_string() const;
	};
	
	/**
	 * A small parser for the Right hand side values of the mexa languages.
	 * This can be
	 * 
	 *    type          example
	 *    ----          -------
	 *    integers      -18     or +15
	 *    doubles       3.1415  or +2.5e8
	 *    booleans      True/False, Yes/No, On/Off  (case insensitive)
	 *    strings       "foo"   or 'foo', or without quotes (see comment below)
	 * 
	 * In any case, single line comments starting with # will be removed.
	 *
	 * This is a lightweight parser which will be quite broad to allow anything.
	 **/
	struct parser {
		symbol leaf;
		std::string input;
		sourcemap src;
		
		parser(symbol leaf, std::string input, sourcemap src) : leaf(leaf),input(input),src(src) {}
		
		static void stripComment(std::string& line);

		/// Parse stuff. Result value is true if output is usable and false if
		/// it could not succeed.
		template<typename T> bool tryCast(T& output) const;
		
		/// An abstract getter for a templated type
		template<typename T> T get_as(std::string type_as_str) const;
		
		/// Read a string at the questioned symbol,
		/// A correctly quoted string.
		std::string get_as_quoted_string() const;
		bool is_a_quoted_string() const;
		
		/// Read a boolean value at the questioned symbol
		bool get_as_bool() const;
		bool is_a_bool() const;
		
		/// Read an integer
		int get_as_int() const;
		bool is_a_int() const;
		
		/// Read a floating point value (double).
		double get_as_double() const;
		bool is_a_double() const;
		
		/**
		 * A parser for regular Mexa-formatted data. This will:
		 *   1) Detect int, double
		 *   2) Strip line comments starting with #
		 *   3) Detect bool's in the mexa variants (True/Yes/On/False/No/Off)
		 *   4) Detect correct strings enclosed in " or '
		 *   5) Fail if nothing is detected with an exception.
		 *      In such a case, the sourcemap is used for error description.
		 **/
		value getValue() const;
		
		/**
		 * A parser for sloppy key value data. This will:
		 *  1) Detect int, double, bool as in fromString()
		 *  2) Understand everything else as string.
		 * Especially, strings shall not be enclosed in " or '.
		 **/
		value getSloppyValue() const;
	};
	
	/**
	* A parser for a simple mexa file. You pass an input stream and we give you
	* a mexafile instance.
	* This parser properly deals with comments and especially requires strings
	* to be enclosed in " or ', as it is defined in the mexa file format.
	* 
	* @param filename_or_desc Used in case of errors to trace back the errnous input.
	* 
	* (One could also come up with a parser for the full mexa file format.)
	**/
	mexafile fromFile(
		std::istream& fh,
		const std::string filename_or_desc="unknown");

	/**
	 * Reads the first line of a (potential) Mexa string and checks for the magic string.
	 * This line is then removed from the istream unless you rewind.
	 **/
	bool hasMagicString(std::istream& fh, const std::string magicString="##mexa");
	
	/**
	* Read a simple-mexa file from an embedded, i.e. encoded, place. Typical formats
	* are base64, base58, base16 as well as `quoted printable' (pass "quotedprintable")
	* where all textual parts are readable and only non-alphanumeric characters are
	* escaped.
	* 
	* This entrypoint is primarily for debugging or in-library use.
	* If you want to read embedded mexa files
	* from ExaHyPE specification files, look at the generic fromSpecfile function.
	**/
	mexafile fromEmbedded(
		const std::string& format,
		const std::string& content,
		const std::string filename_or_desc="unknown");

	/**
	 * Construct a mexafile from a sloppy STL vector with LHS and RHS seperated as
	 * a pair. It will run the parser::getSloppyValue() on the RHS while constructing
	 * the mexafile.
	 * 
	 * After parsing, this function also looks for certain keys in the mexafile which
	 * indicate that there is actually an mexa file embedded in a single string on the
	 * RHS of this ordered map. If so, it will unpack this embedded mexa file.
	 * 
	 * This is the primarily interface to ExaHyPE. We do not allow the
	 * exahype::Parser::ParserView instanc  here for not having a
	 * reference to the ExaHyPE code at this place (seperation of code bases).
	 * 
	 * The usage boils down to as simple as
	 *
	 *  exahype::Parser::ParserView view(parser.getParserView(0)); // obtain it
	 *  mexafile foo = exa::fromSpecfile(view.getAllAsOrderedMap(), "specfile.filename...");
	 *
	 **/
	mexafile fromSpecfile(
		const std::vector<std::pair<std::string, std::string>>& list,
		const std::string source_description="unknown");

} // ns mexa

#endif /* MEXA_META_EXAHYPE_PARFILE_CPP_IMPL_HEADER */
