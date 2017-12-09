#!/usr/bin/python
# Regular Python 2
# A graph approach to mexa.

# This is a clean rewrite/refactor of mexa. 

# dependency
import networkx as nx


# batteries:
import os, re, sys, ast, inspect, argparse, base64, pprint, subprocess
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
tuplize = lambda x: (x,)
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

### end helpers

class term:
	# a unique thing (string). For nxNetwork graph
	# Only one term is equal to all others: The None or "" term
	@classmethod
	def _incrCounter(cls):
		if not hasattr(cls, 'uniq_counter'):
			cls.uniq_counter = 0
		cls.uniq_counter += 1
		return cls.uniq_counter
	def __init__(self,value=None, addr_instance=None):
		"""
		Get a new term. This will be an unique object, i.e. term("foo") != term("foo").
		In order to address an exixting term (mostly for debugging purposes), you can
		write term("foo", 7) if your existing term("foo") has a counter==7.
		"""
		self.value = value
		if addr_instance == None:
			self.counter = self.__class__._incrCounter() if value else 0
		else:
			self.counter = addr_instance
	def isRoot(self):
		return not bool(self.value)
	def copy(self):
		# make a new instance of term which is unrelated to this one.
		return term(self.value)
	def __repr__(self):
		return "term(%s,%d)"%(repr(self.value),self.counter) if self.value else "root"
	def __eq__(self,other):
		#return not self.value and not other.value
		return isa(other,term) and self.counter == other.counter and self.value == other.value
	def __hash__(self):
		return hash((self.value,self.counter))

# define the root symbol for convenience
term.root = term()

class symbol:
	"""
	A symbol is an hierarchical identifier, similar to symbols in LISP and Mathematica.
	Symbols are the atoms of the mexa language.
	Symbols are treated as immutable: There is no method to change them after construction.
	In Mexa, symbols always appear as the LHS of an relation (operation).
	"""

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
	def as_str_tuple(self, include_root=False):
		# Returns the (immutable) tuple of strings holding the elements of the path
		return ("",)+self._path if include_root else self._path
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
		return isa(other,symbol) and self.canonical() == other.canonical()
	def __lt__(self, other): # sortable
		return self.canonical() < other.canonical()


class source:
	"""
	A source is a list of where something comes from.
	Sources just have to yield strings when asked for. That's it.
	"""
	def __init__(self, src=None):
		if isa(src, source):
			self.sources = src.sources
		elif isa(src, list): # shall not catch sourceline instances
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

source.unknown = source()

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
	symbol_alternatives = { '=': 'equals', '<': 'subsets', '<=': 'subsets', '<<': 'include', '+=': 'append' }
	primitives = [int,float,str,unicode,bool]
	opor = "|".join(map(re.escape, symbol_alternatives.keys())) # regex detecting operators
	symb_split = ("::", "/")
	symb_split_or = "|".join(map(re.escape, symb_split))
	symb = r"[a-z_][a-z_0-9:/]*" # regex defining a LHS symbol
	comment = r"#" # comment character
	linecomment = r"(?:" + comment + r".*)?$" # allows a comment until end of line
	stringvarchar = "@" # variable escape character
	stringvarsimple = stringvarchar + "([a-z_0-9]+)" # only alphanumeric
	stringvarcomplex = stringvarchar + "\{("+symb+")\}"
	
	class evalstr:
		def __init__(self, text):
			self.text = text
		def __call__(self, context):
			replmatch = lambda matchobj: context.resolve_symbol(matchobj.group(1))
			text = re.sub(lang.stringvarsimple, replmatch, text, count=0, flags=baseflags)
			text = re.sub(lang.stringvarcomplex, replmatch, text, count=0, flags=baseflags)
			return lang.rhs2python(text, context=context)
		def __repr__(self):
			return "%s(%s)" % (self.__class__.__name__, self.text)

	
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
			# if we already have an operations context, directly evaluate the variables.
			# Otherweise, return the closure.
			promise = lang.evalstr(text)
			return promise(context) if context else promise
		
		# We use python's parser for both understanding the data structure
		# and removing pythonic line comments such as "#".
		try:
			typed = ast.literal_eval(text)
			if isa(typed, lang.primitives):
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
	def rhs2string(cls, rhs):
		if isa(rhs,symbol):
			return rhs.canonical()
		if isa(rhs,lang.evalstr):
			return '"%s"' % rhs.text
		if isa(rhs,lang.primitives):
			return str(rhs)
		else:
			raise ValueError("Bad RHS, cannot transform safely to text: %s" % str(rhs))
	
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
		op = lang.symbol_alternatives[parts.group('op')]

		return Rel(name, value, op, srcline)
	
	@classmethod
	def Rel_to_textline(cls, rel):
		return "%s %s %s" % (rel.l.canonical(), invdict(lang.symbols)[rel.op], cls.rhs2string(rel.r))

class Rel:
	def __init__(self, l, r, op, src=None):
		self.l = l # typically a symbol. May also be a term in some contexts.
		self.r = r # typically something or also a symbol
		self.op = op
		self.src = source(src)
	def __repr__(self):
		return "%s(%s, %s)" % (self.op, self.l, self.r) # omit src for the time being

	def evaluate_rhs(self, *arg):
		if callable(self.r):
			self.r = self.r.value(*arg)
		return self.r
	
	def add_prefix(self, path): # to be updated or removed
		"Return a new operation where lhs is prefixed with path"
		return self.__class__( self.l.add_prefix(path), self.r, self.op, self.src)
	
	def remove_prefix(self, path): # to be updated or removed
		"Return a new operation list where each lhs has a common prefix with path removed"
		return self.__class__( self.l.remove_prefix(path), self.r, self.op, self.src)
	
	# list or tuple idiom
	def __len__(self):
		return 4
	def __iter__(self):
		yield self.l
		yield self.r
		yield self.op
		yield self.src

class mexagraph:
	# the graph has two types of edges
	P = 'path'
	E = 'extends' # or inherits
	
	def __init__(self, edges):
		self.G = nx.DiGraph()
		# add a root
		self.G.add_node(term.root)
		self.from_mexafile(edges)
	def __repr__(self):
		return '%s(%s)' % (self.__class__.__name__, pprint.pformat(self.G.edges(data=True),width=1))
	
	def insert_path(self, path, starting_from=term.root, src=None):
		"""
		Insert edges for a  symbol path, i.e. /a/b/c gets a-[path]->b-[path]->c
		Returns the last edge target which was inserted.
		"""
		if not isa(path,symbol): path = symbol(path)
		if path.isRoot():
			return starting_from
		head, tail = path.base().canonical(), path.except_base()
		child = None
		for existing_child, attr in self.G[starting_from].iteritems():
			if attr['op'] == self.P:
				if existing_child.value == head:
					child = existing_child
					break
		# if did not found the segment, create it:
		if not child:
			child = term(head)
			self.G.add_edge(starting_from, child, op=self.P, src=src)
		return self.insert_path(tail, child, src=src)

	def get_path_down(self, left, starting_from=term(), silent_failure=False): # the original get_path
		"""
		Resolve symbol -> term.
		starting from some term (default: root).
		Resolves down the graph.
		Neccessary because term() instances are unfindable by definition.
		Returns None if nothing found, raises exception if asked for.
		Symbol is always supposed to be rootet at starting_from.
		-> works!
		"""
		if not isa(left,symbol): left = symbol(left)
		if left.isRoot(): # we found our element
			return starting_from
		first = left.base().canonical()
		for right, attr in self.G[starting_from].iteritems():
			if attr['op'] == self.P:
				if right.value == first:
					#print "Looking for %s, Found %s" % (left,first)
					return self.get_path_down(left.except_base(), right)
				else:
					#print "No success: %s != %s" % (right,first)
					pass
		#print "Looked for %s (%s), found nothing" % (left,first)
		if not silent_failure:
			raise ValueError("Symbol %s not found in graph, searching from %s." % (left,starting_from))
		
	def get_path_up(self, symb, right, include_start=True, stack=tuple()):
		"""
		Resolves symbol -> [term], starting from some term.
		# Resolve from a term to a symbol by going back the path
		# symb: Symbol to look for, e.g.  foo/bar
		# right: Term where to start looking at, i.e. biz in bla/boo/biz
		# stack: Internal stack to avoid loops
		# In this example, foo/bar could be found at bla/foo/bar or bla/boo/foo/bar
		# or even bla/boo/biz/foo/bar.
		
		It always returns a list. Elements of this list are tuples.
		First value is target symbol, second value is the stack i.e. path from
		starting symbol to final value (neccessary for understanding double resolutions).
		Ideally, only one list element is returned and the resolution is trivial.
		"""
		
		if not isa(symb,symbol): symb = symbol(symb)
		ret = []
		# test from right itself
		if include_start:
			resolving_term = self.get_path_down(symb, starting_from=right, silent_failure=True)
			if resolving_term:
				ret.append( resolving_term )
		# check in parents of right
		for left, attr in self.G.pred[right].iteritems():
			resolving_term = self.get_path_down(symb, starting_from=left, silent_failure=True)
			if resolving_term:
				ret.append(resolving_term)
			else:
				if left in [rel.l for rel in stack]:
					# run into an equality a-[equals]->b, b-[equals]->a. Do not follow.
					# Equality loops are there by design.
					continue
				else:
					# traverse over the edge
					newstack = stack + tuplize(Rel(left, right, attr['op'], attr['src']))
					ret += self.get_path_up(symb, right=left, stack=newstack)
		return unique_preserve(removeFalse(ret))
	
	# currently UNUSED.
	def resolve_symbol(self, symb, start, silent_failure=False):
		"""
		Returns a symbol to a term.
		"""
		res = self.get_path_up(symb, right=start, include_start=True)
		if len(res) == 1:
			return first(res)
		elif len(res) == 0:
			if silent_failure:
				return None
			else:
				raise ValueError("Symbol %s not found around term %s" % (symb,start))
		else: # len(res)>1
			raise ValueError("Found multiple candidates for %s around term %s: %s" % (symb,start,res))

	def from_mexafile(self, edges):
		"""
		Adds edges to the graph. The rules are:
		  * paths are resolved into the graph and get path edges
		  * equalities and subsets are represented as extends edges
		  * all other operations are not touched and go throught the processing.
		That is, you should make sure if you have operations such as "include" and "append",
		parse them before.
		"""
		
		for l,r,op,src in edges:
			assert isa(l,symbol)
			l = self.insert_path(l, src=src)
			
			if isa(r,symbol):
				target = self.get_path_up(r, l) # check if variable exists
				if len(target) == 1:
					r = first(target) # link to existing
				elif len(target) == 0:
					r = self.insert_path(r, src=src) # insert new node
				else: # len(target) > 1:
					# multiple targets exist. This is bad.
					raise ValueError("Found multiple candidates for %s around term %s: %s" % (r,l,target))
			# => allows full relative variable addressation.
			# => problem: Resolving variables may not yet take future relations into account.
			#    to encounter the problem, the graph would have needed to be setup in a
			#    iterative process with correction steps.
			
			# map equals and subsets together
			if op == 'equals':
				self.G.add_edge(l, r, op=self.E, src=src)
				self.G.add_edge(r, l, op=self.E, src=src) # TODO: Should comment the source.
			elif op == 'subsets':
				self.G.add_edge(l, r, op=self.E, src=src)
			else:
				self.G.add_edge(l, r, op=op, src=src)
	
	# graph to mexa: Without evaluation of any edges
	def to_mexafile(self, left=term.root, left_path=symbol()):
		#assert isa(root, term)
		ret = []
		for right, attr in self.G[left].iteritems():
			#print "%s,%s" % (right,attr)
			if attr['op'] == self.P:
				right_path = symbol(right.value).add_prefix(left_path)# if isa(right,term) else str(right))
				#print "At %s, %s: Visiting %s, %s" % (left, left_path, right, right_path)
				ret += self.to_mexafile(right, right_path)
			else:
				#print "%s: Attr is %s" % (right, attr)
				if isa(right, term):
					right = self.get_path_down(right) # TODO: Resolve the term -> symbol here.
				ret += [ Rel(left_path, right, attr['op'], attr['src']) ]
		return mexafile(ret)
		#return map(unpack(node_to_mexafile), )
		
	# check for undefined symbols
	def check_undef(self, symb=term.root):
		# undefined symbols are defined simply as: A path leaf which has no outoing equalities.
		#[x for x in self.G.nodes_iter() if self.G.out_degree(x)==0 ]
		pass
		# TOD BE DONE.
		
	def evaluate_to_mexafile(self, left=term.root, stack=tuple(), max_rec=10):
		"""
		Correctly resolves all equalities and directed equalities to an oplist aka a mexafile.
		Current Limitations:
		  * Code ignores loops in the equality graph, thus it cannot detect
		    missing definitions. This should be done independently by check_undef.
		"""
		if max_rec == 0:
			raise ValueError("self.Path is too deep. at %s, stack=%s" % (left,stack))
		ret = []
		for right, attr in self.G[left].iteritems():
			if attr['op'] == self.P:
				# put onto the stack:
				newstack = stack + tuplize(Rel(left, right, self.P, attr['src']))
				#right_path = symbol(right.value).add_prefix(left_path)
				#print "At %s, %s: Visiting %s, %s" % (left, left_path, right, right_path)
				ret += self.evaluate_to_mexafile(right, newstack, max_rec=max_rec-1)
			elif attr['op'] == self.E:
				if isa(right, term):
					if right in [rel.l for rel in stack]:
						# run into an equality a-[equals]->b, b-[equals]->a. Do not follow.
						# Equality loops are there by design.
						continue
					else:
						# follow the equals with same path
						newstack = stack + tuplize(Rel(left, right, self.E, attr['src']))
						ret += self.evaluate_to_mexafile(right, newstack, max_rec=max_rec-1)
				else:
					# make an attribute node.
					newstack = stack + tuplize(Rel(left, right, self.E, attr['src']))
					pathlst = [ rel.r.value for rel in newstack if rel.op == self.P ]
					srclst  = [ rel.src     for rel in newstack if rel.op == self.E and not rel.l.isRoot() ]
					#print "At left=%s, newstack=%s I composed path=%s, srclist=%s" % (left,newstack,pathlst,srclst)
					ret += [ Rel(symbol(pathlst), right, "equals", source(srclst)) ] # could also use "define"
			else:
				# TODO: Should instead let all other operations pass.
				raise ValueError("At l=%s,stack=%s, Operation not known: rattr=%s" % (left,stack,attr))
		#if len(ret) == 0:
			# no childs given at this node = undefined!
		#	raise ValueError("Term l=%s,stack=%s lacks a definition\n"%(left,stack))
		return mexafile(ret)
		#return map(unpack(node_to_mexafile), )

class mexafile:
	"""
	Current mexafile
	"""
	
	def __init__(self, ops=None):
		if not ops:
			self.ops = list()
		elif isa(ops,mexafile):
			self.ops = ops.ops # reference
		else:
			self.ops = ops # reference
		# we do not call evaluate() here, do this manually.
		
	def __repr__(self):
		return '%s(%s)' % (self.__class__.__name__, pprint.pformat(self.ops,width=1))
	def __iter__(self):
		return iter(self.ops)
	def __len__(self):
		return len(self.ops)
	
	def add_source(self, something): # chainable
		for op in self.ops:
			op.src.add_source(something)
		return self
	
	def graph(self):
		# There is no mechanism to keep ops and graph in sync.
		if not hasattr(self, '_graph'):
			self._graph = mexagraph(self)
		return self._graph
	
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
	

	def evaluate_symbol(self, symb, evaluator, inplace=True, eliminate=True, max_rec=10):
		"""
		Evaluate a questioned symbol, where symbol is a class, for instance `append`.
		Returns new operations object or does the evaluation inplace.
		"""
		ret_oplist = []
		for rel in self.ops:
			if rel.op == symb:
				new_oplist = evaluator(rel, self)
				ret_oplist += new_oplist
			else:
				ret_oplist.append(rel) # pass throught
		
		if eliminate:
			remaining = mexafile([ rel for rel in ret_oplist if rel.op == symb ])
			if any(remaining):
			# There are still instances of symbol in the oplist and we were
			# asked to eliminate all of them. Recursively call ourselves,
				# expecting that they vanish.
				if max_rec == 0:
					raise ValueError("While trying to evaluate %s, reached maximum number of iterations. The user probably included cyclic links. These symbols remain: %s" % (str(symb),str(remaining)))
				return self.evaluate_symbol(symb, evaluator, inplace=inplace, eliminate=eliminate, max_rec=max_rec-1)
		
		if inplace:
			self.ops = ret_oplist
		else:
			return mexafile(ret_oplist)
		
	def evaluate_all_symbols(self, evaluator, inplace=True):
		ret = evaluator(self)
		if inplace:
			self.ops = ret
		else:
			return mexafile(ret)
	
	def evaluate(self, inplace=True):
		"""
		The evaluation is two-place: First, there is a element-local replacement step.
		Second, there is a global variable resolving step, using the graph.
		"""
		
		def include(op, operations):
			# include a file
			fname = op.evaluate_rhs(self.resolver_for(op))
			if not isa(fname, [str,unicode]):
				raise ValueError("For file inclusion, only strings are supported. In %s" % str(op))
			return mexafile.from_filename(fname).add_prefix(op.l).add_source(op.src)

		def append(op, operations):
			# This algorithm reads as rule: "a/b += c/d" => "a/b/d = c/d"

			if isa(op.r,list):
				# we do support multiple RHS values, i.e. a syntax like
				#  a += b c d  <=> equals(a, sourced([b,c,d], ...))
				return flatten2d([append(o, operations) for o in op.r])
			
			# name of the node to create below the lhs.
			if isa(op.r, symbol):
				name = op.r.node()
			else:
				# come up with some name describing this object
				# TODO: Name should be improved.
				name =  hex(abs(hash( op.r )%2**30))
				
			return mexafile([ Rel(symbol(name).add_prefix(op.l), op.r, 'equals', op.src) ])
		
		def prepare_equal(op, operations):
			if isa(op.r.value, list):
				# lists get created by rhs2python i.e. with something like a = (1,2,3)
				def listToAppends(i,vi): #for i,vi in enumerate(op.r.value):
					if isa(vi,symbol.primitives):
						return Rel(symbol("l%d"%i).add_prefix(op.l), vi, 'equals', op.src.add_source("list expansion"))
					else:
						raise ValueError("self.evaluate_symbol('include', include)Found illegal non-primitive value in RHS of assignment operation "+str(op))
				return map(listToAppends, enumerate(op.r.value))
			else:
				return [op]
		
		def equalities(operations): # caveat, this is global
			# resolves = and <=
			return operations.graph().evaluate_to_mexafile().ops
		
		self.evaluate_symbol('include', include, inplace=inplace)
		self.evaluate_symbol('append', append, inplace=inplace)
		self.evaluate_symbol('equal', prepare_equal, eliminate=False, inplace=inplace)
		self.evaluate_all_symbols(equalities, inplace=inplace)
		# as a last step, should look for inconsistencies or check whether all data
		# have correctly been evaluated.

	def resolve_symbol(self, varname, src=source.unknown):
		"""
		Looks up the value of a variable such as "foo" or "foo/bar" or "foo::bar::baz"
		in the symbols dictionary. In case of errors, src is spilled out.
		"""
		sv = symbol(varname)
		for rel in self.ops:
			if rel.l == sv:
				return rel.r
		raise ValueError("Variable '%s' not defined but used in %s" % (varname,src.verbose()))
	
	def resolver_for(self, rel):
		"Returns a function for resolving a variable. Used for the evalstr() instances."
		return lambda varname: operations.resolve_symbol(rel.r, varname, rel.src)

	def add_prefix(self, path):
		"Return a new operations list where each lhs is prefixed with path"
		if not isa(path,symbol): path = symbol(path)
		return mexafile([ op.add_prefix(path) for op in self.ops ])
	
	def remove_prefix(self, path):
		"Return a new operation list where each lhs has a common prefix with path removed"
		if not isa(path,symbol): path = symbol(path)
		return mexafile([ op.remove_prefix(path) for op in self.ops ])

	def toPlain(self):
		return "\n".join([ lang.Rel_to_textline(rel) for rel in self ])
	
	def toGraph(self, base_filename):
		graph = self.graph()
		dotfilename = base_filename+".dot"
		imgfilename = base_filename+".png"
		nx.nx_agraph.write_dot(graph.G, dotfilename)
		retcode = subprocess.call(["dot", "-Grankdir=LR", "-Tpng", dotfilename, "-o", imgfilename])
		print "Wrote GraphViz output to %s, invoked `dot`, produced %s with exit code %d." % (dotfilename, imgfilename, retcode)


parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
	metavar='FILENAME', help="Input file to read. If no file is given, read from stdin.")
parser.add_argument('--outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
	metavar='FILENAME', help="Output file to write to. If no file is given, write to stdout.")
args = parser.parse_args()

mf = mexafile.from_filehandle(args.infile)
mf.evaluate()

g = mf.graph()

print mf.toPlain()
mf.toGraph("graph")

# render on terminal with
# dot -Grankdir=LR -Tpng graph.dot  -o graph.png

###nx.nx_agraph.write_dot(mf.graph.G, "graph.dot")


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
