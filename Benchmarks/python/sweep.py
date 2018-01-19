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
    result = parseArgument(argv,1) not in ["build","generate"] or \
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

def hashDictionary(dictionary):
    """
    Hash a dictionary.
    """
    chain = ""
    for key,value in dictionary.items():
        chain += key+","+value+";"
    
    result = hashlib.md5(chain.encode()).hexdigest()
    return result
    
def renderSpecificationFile(templateBody,buildParameterDict):
    renderedFile = templateBody
    
    for key,value in buildParameterDict.items():
        renderedFile = renderedFile.replace("{{"+key+"}}", value)
    return renderedFile

if __name__ == "__main__":
    import sys,os
    import configparser
    import subprocess
    import itertools
    import hashlib
    
    if haveToPrintHelpMessage(sys.argv):
        print("sample usage:./sweep.py (build|generate) options.sweep")
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
        templateFileName = workspace["template"]
        exahypeRoot      = workspace["exahype_root"]
        outputPath       = workspace["output_path"]
        
        templateBody = None
        with open(exahypeRoot + "/" + templateFileName, "r") as templateFile:
            templateBody=templateFile.read()
        
        if templateBody!=None:
            if not os.path.exists(exahypeRoot+"/"+outputPath):
                os.makedirs(exahypeRoot+"/"+outputPath)
            if not os.path.exists(exahypeRoot+"/"+outputPath+"/build"):
                os.makedirs(exahypeRoot+"/"+outputPath+"/build")
            
            # build-specific parameters
            if "dimension" not in parameterspace.keys():
                parameterspace["dimension"] = [""]
            if "order" not in parameterspace.keys():
                parameterspace["order"] = [""]
            dimensions = parameterspace["dimension"]
            orders     = parameterspace["order"]
            
            environmentProduct = dictProduct(environmentspace)
            parametersProduct  = dictProduct(parameterspace)
            
            buildParameterDict = list(parametersProduct)[0]
            
            for environmentDict in environmentProduct:
                for key,value in environmentDict.items():
                    os.environ[key]=value
                
                for dimension in dimensions:
                    # make clean
                    print("make clean")
                    process = subprocess.Popen(["make clean"], stdout=subprocess.PIPE, shell=True)
                    (output, err) = process.communicate()
                    process.wait()
                    
                    for order in orders:
                        buildParameterDict["dimension"]=dimension
                        buildParameterDict["order"]    =order
                        
                        buildSpecificationFileBody = renderSpecificationFile(templateBody,buildParameterDict)
                        
                        projectName=workspace["project_name"]
                        buildSpecificationFileName = outputPath + "/build/" + projectName + "-d" + dimension + "-p" + order + ".exahype"
                        
                        with open(exahypeRoot + "/" + buildSpecificationFileName, "w") as buildSpecificationFile:
                            buildSpecificationFile.write(buildSpecificationFileBody)
                        # run toolkit
                        toolkitCommand = "(cd "+exahypeRoot+" && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive "+buildSpecificationFileName+")"
                        print(toolkitCommand);
                        process = subprocess.Popen([toolkitCommand], stdout=subprocess.PIPE, shell=True)
                        (output, err) = process.communicate()
                        process.wait()
                        if "setup build environment ... ok" in str(output):
                            print(toolkitCommand+ " ... ok")
                        else:
                            print("ERROR: "+toolkitCommand + " ... FAILED!",file=sys.stderr)
                            sys.exit()
                        
                        # clean locally & call make
                        gmake_threads=workspace["gmake_threads"]
                        makeCommand="rm -r *.o cipofiles.mk cfiles.mk ffiles.mk kernels; make -j"+gmake_threads
                        print(makeCommand)
                        process = subprocess.Popen([makeCommand], stdout=subprocess.PIPE, shell=True)
                        (output, err) = process.communicate()
                        process.wait()
                        if "build of ExaHyPE successful" in str(output):
                            print(makeCommand+ " ... ok")
                        else:
                            print("ERROR: "+gmakeCommand + " ... FAILED!",file=sys.stderr)
                            sys.exit()
                        
                        projectPath   = workspace["project_path"]
                        oldExecutable = exahypeRoot+"/"+projectPath+"/ExaHyPE-"+projectName
                        newExecutable = exahypeRoot+"/"+outputPath+"/build/ExaHyPE-"+projectName+"-"+hashDictionary(environmentDict)+"-d" + dimension + "-p" + order
                        moveCommand   = "mv "+oldExecutable+" "+newExecutable
                        print(moveCommand)
                        subprocess.call(moveCommand,shell=True)
        else:
            print("ERROR: Couldn\'t open template file: "+workspace["template"])
                    
    elif subprogram == "generate":
        environmentProduct = dictProduct(environmentspace)
        parametersProduct  = dictProduct(parameterspace)
        
        # These hash functions are not robust at all yet.
        # probably have to write my own
        for myTuple in environmentProduct:
          print(hashDictionary(myTuple))
          
        for myTuple in parametersProduct:
          print(hashDictionary(myTuple))
    
