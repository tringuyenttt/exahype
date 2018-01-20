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
    result = parseArgument(argv,1) not in subprograms or \
             parseArgument(argv,2)==None
    for arg in argv:
        result = result or ( arg=="-help" or arg=="-h" )
    return result

def parseEnvironment(config):
    """
    Parse the environment section.
    """
    environmentspace = {}
    if "environment" in config and len(config["environment"].keys()):
        for key, value in config["environment"].items():
            environmentspace[key] = [x.strip() for x in value.split(",")]
        if "SHAREDMEM" not in environmentspace:
            print("ERROR: 'SHAREDMEM' missing in section 'environment'.",file=sys.stderr)
            sys.exit()
    else:
        print("ERROR: Section 'environment' is missing or empty! (Must contain at least 'SHAREDMEM'.)",file=sys.stderr)
        sys.exit()
    
    return environmentspace


def parseParameters(config):
    """
    Parse the parameters section.
    """
    parameterspace = {}
    if "parameters" in config and len(config["parameters"].keys()):
        for key, value in config["parameters"].items():
            parameterspace[key] = [x.strip() for x in value.split(",")]
            
        if "order" not in parameterspace:
            print("ERROR: 'order' missing in section 'parameters'.",file=sys.stderr)
            sys.exit()
        elif "dimension" not in parameterspace:
            print("ERROR: 'dimension' missing in section 'parameters'.",file=sys.stderr)
            sys.exit()
    else:
        print("ERROR: Section 'parameters' is missing or empty! (Must contain at least 'dimension' and 'order'.)",file=sys.stderr)
        sys.exit()
    
    return parameterspace

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
    for key,value in sorted(dictionary.items()):
        chain += key+","+value+";"
    
    result = hashlib.md5(chain.encode()).hexdigest()
    return result
    
def renderSpecificationFile(templateBody,parameterDict,tasks,cores):
    renderedFile = templateBody
    
    for key,value in parameterDict.items():
        renderedFile = renderedFile.replace("{{"+key+"}}", value)
    renderedFile = renderedFile.replace("{{tasks}}", value)
    renderedFile = renderedFile.replace("{{cores}}", value)
    return renderedFile

def build(workspace,machine,environmentspace,parameterspace):
    """
    Build the executables.
    """
    templateFileName = workspace["template"]
    exahypeRoot      = workspace["exahype_root"]
    outputPath       = workspace["output_path"]
    
    templateBody = None
    with open(templateFileName, "r") as templateFile:
        templateBody=templateFile.read()
    
    if templateBody!=None:
        if not os.path.exists(outputPath+"/build"):
            os.makedirs(outputPath+"/build")
        
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
                    
                    buildSpecificationFileBody = renderSpecificationFile(templateBody,buildParameterDict,"1","1")
                    
                    projectName=workspace["project_name"]
                    buildSpecificationFileName = outputPath + "/build/" + projectName + "-d" + dimension + "-p" + order + ".exahype"
                    
                    with open(buildSpecificationFileName, "w") as buildSpecificationFile:
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
                    oldExecutable = projectPath+"/ExaHyPE-"+projectName
                    newExecutable = outputPath+"/build/ExaHyPE-"+projectName+"-"+hashDictionary(environmentDict)+"-d" + dimension + "-p" + order
                    moveCommand   = "mv "+oldExecutable+" "+newExecutable
                    print(moveCommand)
                    subprocess.call(moveCommand,shell=True)
    else:
        print("ERROR: Couldn\'t open template file: "+workspace["template"],file=sys.stderr)

def clean(workspace,subFolder=""):
    outputPath  = workspace["output_path"]
    
    folder = outputPath+"/"+subFolder
    print("rm -r "+folder)
    subprocess.call("rm -r "+folder, shell=True)

def parseCores(jobs,cpus):
    """
    If we encounter "auto" as value, the number of cores is chosen as: 
    total number of cpus (per node) / number of tasks (per node).
    """
    tasks = [x.strip() for x in jobs["tasks"].split(",")]
    cores = [x.strip() for x in jobs["cores"].split(",")]
    if len(cores)==1 and cores[0]=="auto":
        cores = [""]*len(tasks)
        for i,t in enumerate(tasks):
            cores[i] = str(int(int(cpus) / int(t)))
    return cores

def scripts(workspace,machine,environmentspace,parameterspace):
    """
    Generate specification files and job scripts.
    """
    templateFileName = workspace["template"]
    exahypeRoot      = workspace["exahype_root"]
    outputPath       = workspace["output_path"]
    
    jobs       = config["jobs"]
    nodeCounts = [item.strip() for item in jobs["nodes"].split(",")]
    taskCounts = [item.strip() for item in jobs["tasks"].split(",")]
    coreCounts = parseCores(jobs,machine["num_cpus"]);
    
    templateBody = None
    with open(templateFileName, "r") as templateFile:
        templateBody=templateFile.read()
    
    if templateBody!=None:
        subFolder = "scripts"
        
        if not os.path.exists(outputPath+"/"+subFolder):
            os.makedirs(outputPath+"/"+subFolder)
        
        environmentProduct = dictProduct(environmentspace)
        parametersProduct  = dictProduct(parameterspace)
        
        # specification files
        for parametersDict in parametersProduct:
            parametersDictHash = hashDictionary(parametersDict)
            
            for tasks in taskCounts:
                for cores in coreCounts:
                  specificationFileBody = renderSpecificationFile(templateBody,parametersDict,tasks,cores)
                  
                  projectName=workspace["project_name"]
                  specificationFileName = outputPath + "/" + subFolder + "/" + projectName + "-" + parametersDictHash + "-t"+tasks+"-c"+cores+".exahype"
                  
                  with open(specificationFileName, "w") as specificationFile:
                      specificationFile.write(specificationFileBody)
        
        # job scrips
        #for environmentDict in environmentProduct:
        #    environmentDictHash = hashDictionary(environmentDict)



if __name__ == "__main__":
    import sys,os
    import configparser
    import subprocess
    import itertools
    import hashlib
    
    subprograms = ["build","scripts", "cleanBuild", "cleanScripts","clean"]
    
    if haveToPrintHelpMessage(sys.argv):
        print("sample usage:./sweep.py ("+"|".join(subprograms)+") options.sweep")
        sys.exit()
    
    subprogram = parseArgument(sys.argv,1)
    configFile = parseArgument(sys.argv,2)
    
    config = configparser.ConfigParser()
    config.optionxform=str
    config.read(configFile)
    
    workspace = config["workspace"]
    machine   = config["machine"]
    
    environmentspace = parseEnvironment(config)
    parameterspace   = parseParameters(config)
    
    # select subprogram
    if subprogram == "clean":
        clean(workspace)
    elif subprogram == "cleanBuild":
        clean(workspace,"build")
    elif subprogram == "cleanScripts":
        clean(workspace,"scripts")
    elif subprogram == "build":
        build(workspace,machine,environmentspace,parameterspace)
    elif subprogram == "scripts":
        scripts(workspace,machine,environmentspace,parameterspace)
