#!/usr/bin/python
# Regular Python 2
# A graph approach to mexa.

# This is a clean rewrite/refactor of mexa. 

# dependency
import networkx as nx


# batteries:
import os, re, sys, ast, inspect, argparse, base64, pprint
from collections import namedtuple
from argparse import Namespace as namespace
from itertools import izip, islice
#from future.utils import raise_from # no batteries

# helpers and shorthands:

baseflags = re.IGNORECASE
match = lambda pattern,string,flags=0: re.match(pattern,string,baseflags+flags)
unpack = lambda f: lambda p: f(*p)
# The identity function. Don't mix up with pythons internal "id"
idfunc = lambda x: x
# A nicer functional list access thing
first = lambda x: x[0]
last = lambda x: x[-1]
# Removes all None from a list
removeNone = lambda l: [ x for x in l if x is not None ]
removeFalse = lambda l: [ x for x in l if x ]
# removes the list wrapper if it contains only one item
unlistOne = lambda l: l[0] if len(l)==1 else l
# Gives the list of duplicates from a list
dupes = lambda a: [x for n, x in enumerate(a) if x in a[:n]]
# Insert list b at position i in list a, replacing element a[i]
replace_with_list_at = lambda i, a, b: a[:i]+b+a[i+1:]
# Remove common prefix of two strings, i.e. remove_prefix("abcdef","abc")=="def"
remove_comstr_prefix = lambda text, prefix: text[text.startswith(prefix) and len(prefix):]
# flatten a 2d list
flatten2d = lambda l: [item for sublist in l for item in sublist]
# unique items
unique = lambda l: list(set(l))
# unique items while preserve the order
def unique_preserve(seq):
	seen = set()
	seen_add = seen.add
	return [x for x in seq if not (x in seen or seen_add(x))]
# invert dictionary
invdict = lambda d: {v: k for k, v in d.iteritems()}
# an own isninstance function which also allows mapping if types is a list.
# Reads "is a" as the Perl isa function. Mimics isinstance.
def isa(obj, types):
	if isinstance(types, list):
		return any([isa(obj, t) for t in types])
	return isinstance(obj, types)
# Gives an iterator over a dict or list
anykey = lambda obj: obj.keys() if isa(obj,dict) else xrange(len(obj))
# walks a complex dict or list structure
def mapComplex(func, node):
	for key in anykey(node):
		item = node[key]
		if isa(item,dict) or isa(item,list):
			mapComplex(func, item)
		else:
			node[key] = func(item)
	return node
# remove elements from list by indexlist
withoutIndices = lambda lst, indices: [i for j, i in enumerate(lst) if j not in indices]
# raise an exception. Useful from lambdas.
def raise_exception(e):
	raise e
# returns string until first occurance of character, not including the character
untilCharacter = lambda txt, char: first(txt.partition(char))
#
def quoted_printable(s, escape='%'):
	"""
	Returns a quoted printable version of the string s with escape character %, no maximum line length
	(i.e. everything in a single line) and everything non-alphanumeric replaced.
	"""
	# sourcecode inspired by quopri, https://github.com/python/cpython/blob/2.7/Lib/quopri.py
	HEX = '0123456789ABCDEF'
	quote = lambda c: escape + HEX[ord(c)//16] + HEX[ord(c)%16] # quote a single character
	needsquote = lambda c: not ('0' <= c <= '9' or 'a' <= c <= 'z' or 'A' <= c <= 'Z') # Whether character needs to be quoted
	return "".join([ quote(c) if needsquote(c) else c for c in s])
#
def window(seq, n=2):
	"Returns a sliding window (of width n) over data from the iterable"
	"   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
	it = iter(seq)
	result = tuple(islice(it, n))
	if len(result) == n:
		yield result
	for elem in it:
		result = result[1:] + (elem,)
		yield result
#
class NamedFunction:
	"Readable names for functions"
	def __init__(self, name, f):
		self.f = f
		self.name = name
	def __call__(self, *args, **kwargs):
		return self.f(*args, **kwargs)
	def __repr__(self):
		return self.name


class term:
	# a unique thing (string). For nxNetwork graph
	# Only one term is equal to all others: The None or "" term
	@classmethod
	def _incrCounter(cls):
		if not hasattr(cls, 'uniq_counter'):
			cls.uniq_counter = 0
		cls.uniq_counter += 1
		return cls.uniq_counter
	def __init__(self,value=None):
		self.value = value
		self.counter = self.__class__._incrCounter() if value else 0
	def copy(self):
		# make a new instance of term which is unrelated to this one.
		return term(self.value)
	def __repr__(self):
		return "term%s(%s)"%(self.counter,repr(self.value)) if self.value else "root"
	def __eq__(self,other):
		#return not self.value and not other.value
		return self.counter == other.counter and self.value == other.value
	def __hash__(self):
		return hash((self.value,self._incrCounter))

class symbol:
	"""
	A symbol is an hierarchical identifier, similar to symbols in LISP and Mathematica.
	Symbols are the atoms of the mexa language.
	Symbols are treated as immutable: There is no method to change them after construction.
	In Mexa, symbols always appear as the LHS of an relation (operation).
	"""
	
	# primitive types which we allow to store on the RHS.
	primitives = [int,float,str,unicode,bool]
	
	def __init__(self,name=''):
		""""
		Creates a symbol from a name such as Foo/Bar/Baz or a path such
		as ['Foo','Bar','Baz']. Without argument, this gives the root symbol
		symbol().
		Since a symbol is immutable, we can store a string and list version at the
		same time for lookup efficiency.
		"""
		if name == None: name = '' # whatever
		if isa(name,(list,tuple)):
			self._path = map(str.lower, map(str.strip, name))
		elif isa(name, symbol):
			self._path = name._path # kind of copy constructor, but only references to other list.
		else:
			name = str(name).lower() # in order to get case-insensitive
			self._path = map(str.strip, re.split(lang.symb_split_or, name))
		# in case of roots and merged paths and so on
		self._path = removeFalse(self._path)
		# make the path immutable
		self._path = tuple(self._path)
		# Compute a string representation, useful for path computations
		self._canonical = "/".join(self._path)
		
	def canonical(self):
		"Canonical string representation with a specific seperator"
		return self._canonical
	def __repr__(self):
		return "symbol(%s)" % self.canonical()
	def isRoot(self):
		return len(self._path) == 0
	
	def base(self): # first non-root element of path
		return symbol(self._path[0]) if not self.isRoot() else symbol()
	def except_base(self): # everything below the root, except the first non-root element
		return symbol(self._path[1:]) if len(self._path)>1 else symbol()
	def node(self): # last element of path
		if self.isRoot():
			return self
		return symbol(self._path[-1])
	def parent(self):
		if self.isRoot():
			raise ValueError("Symbol is already root")
		return symbol(self._path[:-1])
	def ancestors(self):
		"""
		Lists all parenting symbols, for instance for /a/b/c it is the
		list [/, /a, /a/b].
		"""
		return [symbol(self._path[:i]) for i,_ in enumerate(self._path)]
	def parts_as_term(self):
		"""
		Like ancestors, just weird. For instance /a/b/c it is
		[root, term('a'), term('b'), term('c')]
		"""
		return [term()] + [term(p) for p in self._path]
	def node_as_term(self):
		if self.isRoot():
			return term()
		return term(self._path[-1])
	def base_as_term(self):
		return term(self._path[0])
	
	def add_prefix(self, other, inplace=False): # should be named add_prefix
		"Returns a new symbol wich is prefixed by self."
		if not isa(other,symbol): other = symbol(other)
		newpath = other._path + self._path;
		if inplace:
			self._path = newpath
			return self
		else:	return symbol(newpath)
	
	def remove_prefix(self, other, inplace=False):
		"Returns a new symbol which has removed the common prefix"
		if not isa(other,symbol): other = symbol(other)
		newpath = remove_comstr_prefix(self.canonical(), other.canonical())
		if inplace:
			self._path = newpath
			return self
		else:	return symbol(newpath)

	# allow symbols to be dict keys:
	def __hash__(self):
		return hash(tuple(self._path))
	def __eq__(self, other):
		return self.canonical() == other.canonical()
	def __lt__(self, other): # sortable
		return self.canonical() < other.canonical()

### end helpers

class source:
	"""
	A source is a list of where something comes from.
	Sources just have to yield strings when asked for. That's it.
	"""
	def __init__(self, src=None):
		if isa(src, source):
			self.sources = src
		elif src==None:
			self.sources = ["unknown"]
		else:
			self.sources = [src]
	def __repr__(self): # todo: make nicer
		return "%s(%s)" % (self.__class__.__name__,self.sources_as_str())
	def add_source(self, src):
		self.sources += [src]
	def sources_as_str(self):
		if len(self.sources)==1:
			return self.sources[0]
		else:
			# TODO: Make this nicer
			return " -> ".join(map(str, self.sources))


class sourceline(namedtuple('sourceline', 'fname linenum text')):
	"""
	Sourceline represents the source of something. It is a key class to allow
	transparent traceback of operations back to their source.
	"""
	
	def __repr__(self):
		return '%s:%s' % (self.fname, self.linenum)
	def verbose(self):
		return '%s line %s: `%s`' % (self.fname, self.linenum, self.text.strip())

	@classmethod
	def from_python(cls):
		"Create a sourceline instance from a python caller position."
		previous_frame = inspect.currentframe().f_back
		(filename, line_number, function_name, lines, index) = inspect.getframeinfo(previous_frame)
		return cls(filename, line_number, lines[0])
	
	class sourcemap:
		"""
		Allows to map an iterable with a given source name to inject sourceline objects.
		"""
		def __init__(self,iterable,source_name=None):
			self.iterable = iterable
			# todo: test with StringIO
			if not source_name:
				try:
					self.source_name = iterable.name
				except AttributeError:
					self.source_name = str(iterable) # hopefully short description
			else:
				self.source_name = source_name
		def iter(self):
			"Returns an iterator over the file"
			offset = 1 # python starts with 0 for line counting but humans start with 1
			return (sourceline(self.source_name,linenum+offset,text) for linenum, text in enumerate(self.iterable))
		def map(self,func):
			return map(func, self.iter())


# Regexps for defining the language
class lang:
	symbols = { '=': 'equals', '<=': "subsets", '<<': "include", '+=': "append" }
	opor = "|".join(map(re.escape, symbols.keys())) # regex detecting operators
	symb_split = ("::", "/")
	symb_split_or = "|".join(map(re.escape, symb_split))
	symb = r"[a-z_][a-z_0-9:/]*" # regex defining a LHS symbol
	comment = r"#" # comment character
	linecomment = r"(?:" + comment + r".*)?$" # allows a comment until end of line
	stringvarchar = "@" # variable escape character
	stringvarsimple = stringvarchar + "([a-z_0-9]+)" # only alphanumeric
	stringvarcomplex = stringvarchar + "\{("+symb+")\}"
	
	@classmethod
	def rhs2python(cls, text, context=None):
		"""
		Parses a right hand side into int, float, boolean, symbol, whatever.
		This does return functions in some cases. They are supposed to be executed with
		an operations instance then for variable substitution.
		If this returns a list or if the returned and subsequently evaluated function
		returns a list, it is subject to the caller to deal with this list accordingly.
		
		context shall be a VariableContext instance or None.
		"""
		text = text.strip()
		
		# as a kind of preprocessor, support "@foo"-like variables. In case they
		# are detected (we uglily try to detect them not in comments), we return
		# instead a promise function for later evaluation with a correct variable
		# context.
		text_decommented = first(text.partition(lang.comment))
		if re.search(lang.stringvarsimple, text_decommented) or re.search(lang.stringvarcomplex, text_decommented):
			#import ipdb; ipdb.set_trace()
			def promise_rhs2python(context, text=text): # text=text is a workaround for py2 missing closures
				replmatch = lambda matchobj: context.resolve_symbol(matchobj.group(1))
				text = re.sub(lang.stringvarsimple, replmatch, text, count=0, flags=baseflags)
				text = re.sub(lang.stringvarcomplex, replmatch, text, count=0, flags=baseflags)
				return rhs2python(text, context=context)
			# if we already have an operations context, directly evaluate the variables.
			# Otherweise, return the closure.
			return promise_rhs2python(context) if context else NamedFunction("evalstr(%s)" % text, promise_rhs2python)
		
		# We use python's parser for both understanding the data structure
		# and removing pythonic line comments such as "#".
		try:
			typed = ast.literal_eval(text)
			if isa(typed,symbol.primitives):
				return typed
			elif isa(typed,[tuple,list]):
				# We said we support lists.
				return list(typed)
			else:
				# we probably have a created a weird type, for instance we
				# parsed a dictionary, which we do not want to support.
				# Instead, thus, return a text. This may not correctly strip comments,
				# thought.
				return text
		except (ValueError, SyntaxError):
			# This is most likely a string or something weird
			
			if re.match(r"^(Yes|True|On)"+lang.linecomment, text, re.IGNORECASE):
				return True
			if re.match(r"^(No|False|Off)"+lang.linecomment, text, re.IGNORECASE):
				return False
			
			# A symbol or a list of symbols.
			symbs = re.match(r"^("+lang.symb+")(\s+"+lang.symb+")*"+lang.linecomment, text, re.IGNORECASE)
			# If only one symbol is detected, return the symbol. Otherwise
			# return the list of symbols which has to be understood correctly.
			if symbs:
				return unlistOne(map(symbol, removeNone(symbs.groups())))
			
			return text
		return text
	
	@classmethod
	def Rel_from_textline(cls, srcline): # was textline2edge
		"""
		Parses a single line. Input is a sourceline object. You can test lines quickly with
		> line2operation(sourceline.from_unknown("foo=10"))
		"""
		if re.match(r"^"+lang.comment+"|^\s*$", srcline.text): # comments
			return None

		parts = re.match(r"^(?P<lhs>(?:"+lang.symb+"|))\s*(?P<op>(?:"+lang.opor+"))\s*(?P<rhs>.+)\s*$", srcline.text, re.IGNORECASE)
		if not parts:
			raise ValueError("Don't understand line %d in %s: '%s'"%(srcline.linenum,srcline.fname,srcline.text))

		# split lhs into name path
		name = symbol(parts.group('lhs'))
		# parse rhs as symbol or whatever
		value = lang.rhs2python(parts.group('rhs'))
		op = lang.symbols[parts.group('op')]

		return Rel(name, value, op, srcline)

class Rel:
	def __init__(self, l, r, op, src=None):
		assert isa(l,symbol)
		self.l = l # typically a symbol
		self.r = r # typically something or also a symbol
		self.op = op
		self.src = source(src)
	def __repr__(self):
		return "%s(%s, %s)" % (self.op, self.l, self.r) # omit src for the time being

	def evaluate_rhs(self, *arg):
		if callable(self.r):
			self.r = self.r.value(*arg)
		return self.r
	
	def add_prefix(self, path):
		"Return a new operation where lhs is prefixed with path"
		return self.__class__( self.l.add_prefix(path), self.rhs)
	
	def remove_prefix(self, path):
		"Return a new operation list where each lhs has a common prefix with path removed"
		return self.__class__( self.l.remove_prefix(path), self.rhs)
	
	# list or tuple idiom
	def __len__(self):
		return 4
	def __iter__(self):
		yield self.l
		yield self.r
		yield self.op
		yield self.src

class mexagraph:
	def __init__(self, edges):
		self.G = nx.DiGraph()
		self.from_mexafile(edges)
	
	def insert(self, path):
		"Returns the last edge target which was inserted"
		assert isa(path,symbol)
		ret = term()
		for a,b in window(path.parts_as_term()):
			self.G.add_edge(a, b, op='path')
			ret = b
		return ret

	def get_term(self, left, starting_from=term()):
		# Resolve from a symbol path, starting from root, to a term
		# Neccessary because term() instances are unfindable by definition
		# returns None if not findable
		if not isa(left,symbol): left = symbol(left)
		if left.isRoot(): # we found our element
			return starting_from
		first = left.base().canonical()
		for right, attr in self.G[starting_from].iteritems():
			if attr['op'] == 'path':
				if right.value == first:
					print "Looking for %s, Found %s" % (left,first)
					return self.get_term(left.except_base(), right)
				else:
					print "No success: %s != %s" % (right,first)
		print "Looked for %s (%s), found nothing" % (left,first)
		
	def get_path(self, right):
		# Resolve from a term to a symbol
		return right # TODO
	
	def from_mexafile(self, edges):
		for l,r,op,src in edges:
			assert isa(l,symbol)
			l = self.insert(l) # include the path
			if isa(r,symbol): r = self.insert(r)
			# include the actual equality
			self.G.add_edge(l, r, op=op, src=src)
	
	# graph to mexa, again
	##### TODO:
	##### At this step, evaluate the graph (is requested).
	##### Then, we do everything inplace: += << = <=
	#####
	def to_mexafile(self, left=term(), left_path=symbol()):
		#assert isa(root, term)
		ret = []
		for right, attr in self.G[left].iteritems():
			#print "%s,%s" % (right,attr)
			if attr['op'] == 'path':
				right_path = symbol(right.value).add_prefix(left_path)# if isa(right,term) else str(right))
				#print "At %s, %s: Visiting %s, %s" % (left, left_path, right, right_path)
				ret += self.to_mexafile(right, right_path)
			else:
				#print "%s: Attr is %s" % (right, attr)
				if isa(right, term):
					right = self.get_path(right)
				ret += [ Rel(left_path, right, attr['op'], attr['src']) ]
		return ret
		#return map(unpack(node_to_mexafile), )
	
	# evaluating is easier
	def evaluate(self):
		for u, v, optype in self.G.edges(data='type'):
			if optype == '+=':
				# addition abbreviation, simple.
				if isa(v,term):
					intermediate = v.copy()
				else:
					# v is a primitive. TODO: Come up with a good name.
					# To do so, get number of total children or so.
					intermediate_name = term("todo")
				self.G.add_edge(u, intermediate, type='path')
				self.G.add_edge(intermediate, v, type='=')
				# TODO: remove the current edge
				
			# TODO: Also loop correctly over all optypes
			elif optype == '=':
				#self.G.add_edge(u, 
				pass # TODO
				

class mexafile:
	"""
	Current mexafile
	"""
	
	def __init__(self, edges=[]):
		# store edges only for debugging
		self.edges = edges # this is by intention not a copy, but just a reference
		# this is where we want to work on
		#self.G = edges2graph(edges)
		
	def __repr__(self):
		return '%s(%s)' % (self.__class__.__name__, pprint.pformat(self.edges,width=1))
	
	@classmethod
	def from_filehandle(cls, fh):
		# The ordered list of operations as they appear in the file
		return cls(removeNone(sourceline.sourcemap(fh).map(lang.Rel_from_textline)))
		# Prepend the list with the rootbase
		# self.ops = self.rootbase + self.ops # no more.

	@classmethod
	def from_filename(cls, fname):
		"Create an instance from a filename instead of filehandle"
		with open(fname, 'r') as fh:
			return cls.from_filehandle(fh)
	
	@classmethod
	def from_textlines(cls, oplines, source_name="unoriginated textlines"):
		"""
		* List of single lines: Are parsed, source_name used as source name.
		oplines must be iterable (i.e. list of strings or an open file handle).
		"""
		if not isa(oplines, list):
			oplines = [oplines]
		oplist = removeNone(sourceline.sourcemap(oplines, source_name).map(lang.Rel_from_textline))
		return cls(oplist)
	
	def evaluate(self):
		#### this is superseded. We do *everything* inline when evaluating.
		def include(op, operations):
			# include a file
			fname = op.evaluate_rhs(VariableContext(op,operations))
			if not isa(fname, [str,unicode]):
				raise ValueError("For file inclusion, only strings are supported. In %s" % str(op))
			return operations.from_filename(fname).add_prefix(op.l)

		def append(op, operations):
			# This algorithm reads as rule: "a/b += c/d" => "a/b/d = c/d"
			
			if isa(op.r,list):
				# we do support multiple RHS values, i.e. a syntax like
				#  a += b c d  <=> equals(a, sourced([b,c,d], ...))
				return flatten2d([append(o, operations) for o in op.r])
			
			# name of the node to create below the lhs. This is by
			# definition the node name.
			if not isa(op.r,symbol):
				raise ValueError("Only support symbol for appending, got %s in %s" % (type(op.r), op))
			# todo: could also support appending
			name = op.r.node()
			return operations([ equals(op.lhs.add_prefix(name), op.rhs) ])



parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
	metavar='FILENAME', help="Input file to read. If no file is given, read from stdin.")
parser.add_argument('--outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
	metavar='FILENAME', help="Output file to write to. If no file is given, write to stdout.")
args = parser.parse_args()

mf = mexafile.from_filehandle(args.infile)

mG=mexagraph(mf.edges)
print mG.to_mexafile()

# render on terminal with
# dot -Grankdir=LR -Tpng graph.dot  -o graph.png
nx.nx_agraph.write_dot(mG.G, "graph.dot")


#append_oplist = removeNone(args.expressions)
#if len(append_oplist):
#	print "Parsing: " + str(append_oplist)
#	mf.append(append_oplist, "Command line argument expression")

# could do here also --env which then adds all or requested environment variables below env/

# appending stuff to any file, as a root
#mf.ops += [
	# mexa = path to the current script directory
	# Rel(symbol("mexa"), os.path.dirname(os.path.realpath(__file__), "equals", sourceline.from_python())),
#]

#sys.exit() # stop here for the time being

# todo: Could offer here also an option to inject various or all environment
# variables, just by passing "--env". Could setup manual assign() operations
# similar to mexafile.rootbase.
# Could also get rid of "let" statements by moving the rootbase outside the class.

#mf.evaluate()
