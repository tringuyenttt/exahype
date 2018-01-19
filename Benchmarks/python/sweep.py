#!/usr/bin/python3
"""
.. module:: sweep
  :platform: Unix, Windows, Mac
  :synopsis: Generate benchmark suites for ExaHyPE.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>, 

:synopsis: Generate benchmark suites for ExaHyPE.
"""
def parseArgument(argv,i):
    if i<len(argv):
        return argv[i]
    else:
        return None

def haveToPrintHelpMessage(argv):
    """
    Check if we have to print a help message.
    """
    result = parseArgument(argv,1)==None or \
             parseArgument(argv,2)==None
    for arg in argv:
        result = result or ( arg=="-help" or arg=="-h" )
    return result

def parseCores(jobs,cpus):
    """
    If we encounter "auto" as value, the number of cores is chosen as: 
    total number of cpus (per node) / number of tasks (per node).
    """
    cores = jobs["cores"].split(",");
    if len(cores)==1 and cores[0]=="auto":
        cores = [""]*len(tasks)
        for i,t in enumerate(tasks):
            cores[i] = int(int(cpus) / int(t))
    return cores

def dictProduct(dicts):
    """
    Computes the Cartesian product of a dictionary of lists as 
    a list of dictionaries.
    
    Gladly copied this code from:
    https://stackoverflow.com/questions/5228158/cartesian-product-of-a-dictionary-of-lists
    
    Example input:
    options = {"number": [1,2,3], "color": ["orange","blue"] }
    
    Example output:
    [ {"number": 1, "color": "orange"},
      {"number": 1, "color": "blue"},
      {"number": 2, "color": "orange"},
      {"number": 2, "color": "blue"},
      {"number": 3, "color": "orange"},
      {"number": 3, "color": "blue"}
    ]
    """
    return (dict(zip(dicts, x)) for x in itertools.product(*dicts.values()))

usedHashes = []

def hashDictionary(dictionary):
    """
    Hash a dictionary
    """
    chain = ""
    for key,value in dictionary.items():
        chain += key+","+value+";"
    
    result = hashlib.md5(chain.encode()).hexdigest()
    while result in usedHashes:
        result = hashlib.md5(chain.encode()).hexdigest()
    usedHashes.append(result)
    return result

if __name__ == "__main__":
    import sys,os
    import configparser
    from subprocess import call
    import itertools
    import hashlib
    
    if haveToPrintHelpMessage(sys.argv):
        print("sample usage: python3 sweep (setup|build|generate) options.sweep")
        print("hint: you might want to hide python3 sweep behind an alias, e.g. alias sweep=\"python <mypath>/sweep.py\"")
        sys.exit()
    
    subprogram = parseArgument(sys.argv,1)
    configFile = parseArgument(sys.argv,2)
    
    config = configparser.ConfigParser()
    config.optionxform=str
    config.read(configFile)
    
    workspace = config["workspace"]
    machine   = config["machine"]
    
    jobs  = config["jobs"]
    nodes = jobs["nodes"].split(",");
    tasks = jobs["tasks"].split(",");
    cores = parseCores(jobs,machine["num_cpus"]);
    #print(nodes); print(tasks); print(cores)
    
    # environment
    environmentspace = {}
    if "environment" in config and len(config["environment"].keys()):
        for key, value in config["environment"].items():
            environmentspace[key] = value.split(",")
    else:
        environmentspace["DUMMY_VAR"] = [""] # We will later on have a loop nest; we thus need at least one element
    
    # parameters
    parameterspace = {}
    if "parameters" in config and len(config["parameters"].keys()):
        for key, value in config["parameters"].items():
            parameterspace[key] = value.split(",")
    else:
        parameterspace["DUMMY_VAR"] = [""]
    
    # select subprogram
    if subprogram == "build":
        # build-specific parameters
        if "dimension" not in parameterspace.keys():
            parameterspace["dimension"] = ["-"]
        if "order" not in parameterspace.keys():
            parameterspace["order"] = ["-"]
        
        dimensions = parameterspace["dimension"]
        orders     = parameterspace["order"]
        
        environmentProduct = dictProduct(environmentspace)
        parametersProduct  = dictProduct(parameterspace)
        
        firstParameterTuple = list(parametersProduct)[0]
        # put parameters in spec
        # create compilable instance of the specfile template
        
        for myTuple in environmentProduct:
            for key,value in myTuple.items():
                os.environ[key]=value
                
            for d in dimensions:
                for p in orders:
                    print(str(myTuple.items())+","+d+","+p)
                    # call("make ",shell=True)
    elif subprogram == "generate":
        environmentProduct = dictProduct(environmentspace)
        parametersProduct  = dictProduct(parameterspace)
        
        # These hash functions are not robust at all yet.
        # probably have to write my own
        for myTuple in environmentProduct:
          print(hashDictionary(myTuple))
          
        for myTuple in parametersProduct:
          print(hashDictionary(myTuple))
    
