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
    result = parseArgument(argv,2) not in subprograms or \
             parseArgument(argv,1)==None
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

def clean(general,subFolder=""):
    exahypeRoot = general["exahype_root"]
    outputPath  = general["output_path"]
    
    folder = exahypeRoot+"/"+outputPath+"/"+subFolder
    print("rm -r "+folder)
    subprocess.call("rm -r "+folder, shell=True)

def renderSpecificationFile(templateBody,parameterDict,tasks,cores):
    renderedFile = templateBody
    
    context = dict(parameterDict)
    context["tasks"] = tasks
    context["cores"] = cores
    
    consistent = True
    for key in parameterDict:
        if "{{"+key+"}}" not in templateBody:
            consistent = False
            print("ERROR: parameter '{{"+key+"}}' not found in spec file template!",file=sys.stderr)
    if not consistent:
        sys.exit()
    
    for key,value in parameterDict.items():
        renderedFile = renderedFile.replace("{{"+key+"}}", value)
    
    return renderedFile

def build(general,environmentspace,parameterspace):
    """
    Build the executables.
    """
    templateFileName = general["spec_template"]
    exahypeRoot      = general["exahype_root"]
    outputPath       = general["output_path"]
    subFolder        = "build"
    
    templateBody = None
    with open(exahypeRoot+"/"+templateFileName, "r") as templateFile:
        templateBody=templateFile.read()
    
    if templateBody!=None:
        if not os.path.exists(exahypeRoot+"/"+outputPath+"/"+subFolder):
            print("create directory:"+exahypeRoot+"/"+outputPath+"/"+subFolder)
            os.makedirs(exahypeRoot+"/"+outputPath+"/"+subFolder)
        
        dimensions = parameterspace["dimension"]
        orders     = parameterspace["order"]
        buildParameterDict = list(dictProduct(parameterspace))[0]
        
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
                    
                    projectName=general["project_name"]
                    buildSpecFileName = outputPath+"/"+subFolder+"/"+projectName+"-d"+dimension+"-p"+order+".exahype"
                    
                    with open(buildSpecFileName, "w") as buildSpecificationFile:
                        buildSpecificationFile.write(buildSpecificationFileBody)
                    # run toolkit
                    toolkitCommand = "(cd "+exahypeRoot+" && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive "+buildSpecFileName+")"
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
                    make_threads=general["make_threads"]
                    makeCommand="rm -r *.o cipofiles.mk cfiles.mk ffiles.mk kernels; make -j"+make_threads
                    print(makeCommand)
                    process = subprocess.Popen([makeCommand], stdout=subprocess.PIPE, shell=True)
                    (output, err) = process.communicate()
                    process.wait()
                    if "build of ExaHyPE successful" in str(output):
                        print(makeCommand+ " ... ok")
                    else:
                        print("ERROR: "+makeCommand + " ... FAILED!",file=sys.stderr)
                        sys.exit()
                    
                    oldExecutable = exahypeRoot + "/" + projectPath+"/ExaHyPE-"+projectName
                    newExecutable = exahypeRoot + "/" + outputPath+"/"+subFolder + "/ExaHyPE-"+projectName+"-"+hashDictionary(environmentDict)+"-d" + dimension + "-p" + order
                    moveCommand   = "mv "+oldExecutable+" "+newExecutable
                    print(moveCommand)
                    subprocess.call(moveCommand,shell=True)
    else:
        print("ERROR: Couldn\'t open template file: "+templateFileName,file=sys.stderr)

def parseCores(jobs):
    """
    If we encounter "auto" as value, the number of cores is chosen as: 
    total number of cpus (per node) / number of tasks (per node).
    """
    cpus = jobs["num_cpus"]
    tasks = [x.strip() for x in jobs["tasks"].split(",")]
    cores = [x.strip() for x in jobs["cores"].split(",")]
    if len(cores)==1 and cores[0]=="auto":
        cores = [""]*len(tasks)
        for i,t in enumerate(tasks):
            cores[i] = str(int(int(cpus) / int(t)))
    return cores

def renderJobScript(templateBody,environmentDict,parameterDict,jobs,
                    jobName,jobFileName,outputFileName,errorFileName,appName,specFileName,
                    nodes,tasks,cores,run):
    """
    Render a job script.
    """
    renderedFile = templateBody
    
    context = {}
    # mandatory
    context["nodes"] = nodes
    context["tasks"] = tasks
    context["cores"] = cores
    context["job_name"]    = jobName
    context["output_file"] = outputFileName
    context["error_file"]  = errorFileName
    
    context["environment"] = json.dumps(environmentDict).replace("\"","\\\"")
    context["parameters"]  = json.dumps(parameterDict).replace("\"","\\\"")
    
    context["job_file"]    = jobFileName
    context["app"]        = appName
    context["spec_file"]  = specFileName
    
    consistent = True
    for key in context:
        if "{{"+key+"}}" not in templateBody:
            consistent = False
            print("ERROR: parameter '{{"+key+"}}' not found in job script template!",file=sys.stderr)
    if not consistent:
        sys.exit()
    
    # optional
    context["mail"]  = jobs["mail"]
    context["time"]  = jobs["time"]
    context["ranks"] = str(int(nodes)*int(tasks))
    
    for key,value in context.items():
        renderedFile = renderedFile.replace("{{"+key+"}}", value)
    
    return renderedFile

def generateScripts(general,environmentspace,parameterspace):
    """
    Generate spec files and job scripts.
    """
    exahypeRoot          = general["exahype_root"]
    outputPath           = general["output_path"]
    buildFolder          = "build"
    scriptsFolder        = "scripts"
    resultsFolder        = "results"
    
    jobs       = config["jobs"]
    nodeCounts = [x.strip() for x in jobs["nodes"].split(",")]
    taskCounts = [x.strip() for x in jobs["tasks"].split(",")]
    coreCounts = parseCores(jobs);
    runs       = int(jobs["runs"])
    
    specFileTemplatePath  = general["spec_template"]
    jobScriptTemplatePath = general["job_template"]
    
    specFileTemplate  = None
    with open(exahypeRoot+"/"+specFileTemplatePath, "r") as templateFile:
        specFileTemplate=templateFile.read()
    if specFileTemplate is None:
        print("ERROR: Couldn\'t open template file: "+specFileTemplatePath,file=sys.stderr)
        sys.exit()
        
    jobScriptTemplate = None
    with open(exahypeRoot+"/"+jobScriptTemplatePath, "r") as templateFile:
        jobScriptTemplate=templateFile.read()
    if jobScriptTemplate is None:
        print("ERROR: Couldn\'t open template file: "+jobScriptTemplatePath,file=sys.stderr)
        sys.exit()
        
    if not os.path.exists(exahypeRoot+"/"+outputPath+"/"+scriptsFolder):
        os.makedirs(exahypeRoot+"/"+outputPath+"/"+scriptsFolder)
    if not os.path.exists(exahypeRoot+"/"+outputPath+"/"+resultsFolder):
        os.makedirs(exahypeRoot+"/"+outputPath+"/"+resultsFolder)
    
    # spec files
    for parameterDict in dictProduct(parameterspace):
        parameterDictHash = hashDictionary(parameterDict)
        
        for tasks in taskCounts:
            for cores in coreCounts:
              specFileBody = renderSpecificationFile(specFileTemplate,parameterDict,tasks,cores)
              
              projectName=general["project_name"]
              specFileName = exahypeRoot + "/" + outputPath + "/" + scriptsFolder + "/" + projectName + "-" + parameterDictHash + "-t"+tasks+"-c"+cores+".exahype"
              
              with open(specFileName, "w") as specFile:
                  specFile.write(specFileBody)
    
    # check if required builts exist
    for environmentDict in dictProduct(environmentspace):
        environmentDictHash = hashDictionary(environmentDict)
        for parameterDict in dictProduct(parameterspace):
            dimension = parameterDict["dimension"]
            order     = parameterDict["order"]
            appName   = exahypeRoot + "/" + outputPath+"/"+buildFolder + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-d" + dimension + "-p" + order
            
            if not os.path.exists(appName):
                print("WARNING: required application '"+appName+"' does not exist!",file=sys.stderr)
    
    # generate job scrips
    for run in range(0,runs):
        for nodes in nodeCounts:
            for tasks in taskCounts:
                for cores in coreCounts:
                    for environmentDict in dictProduct(environmentspace):
                        environmentDictHash = hashDictionary(environmentDict)
                        
                        for parameterDict in dictProduct(parameterspace):
                            parameterDictHash = hashDictionary(parameterDict)
                            
                            dimension = parameterDict["dimension"]
                            order     = parameterDict["order"]
                            appName   = exahypeRoot + "/" + outputPath+"/"+buildFolder + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-d" + dimension + "-p" + order
                            
                            
                            specFileName = exahypeRoot + "/" + outputPath + "/" + scriptsFolder + "/" + projectName + "-" + \
                                           parameterDictHash + "-t"+tasks+"-c"+cores+".exahype"
                            
                            jobName        = projectName + "-" + environmentDictHash + "-" + parameterDictHash + \
                                             "-n" + nodes + "-t"+tasks+"-c"+cores+"-r"+str(run)
                            jobFilePrefix  = exahypeRoot + "/" + outputPath + "/" + scriptsFolder + "/" + jobName
                            jobFileName    = jobFilePrefix + ".job"
                            outputFileName = exahypeRoot + "/" + outputPath + "/" + resultsFolder + "/" + jobName + ".out"
                            errorFileName  = exahypeRoot + "/" + outputPath + "/" + resultsFolder + "/" + jobName + ".err"
                            
                            jobScriptBody = renderJobScript(jobScriptTemplate,environmentDict,parameterDict,jobs,
                                                            jobName,jobFileName,outputFileName,errorFileName,appName,specFileName,
                                                            nodes,tasks,cores,run)
                            with open(jobFileName, "w") as jobFile:
                                jobFile.write(jobScriptBody)

def submitJobs(general):
    exahypeRoot          = general["exahype_root"]
    outputPath           = general["output_path"]
    scriptsFolder        = "scripts"
    
    jobScriptFolder = exahypeRoot+"/"+outputPath+"/"+scriptsFolder
    if not os.path.exists(jobScriptFolder):
        print("ERROR: Job script folder '"+jobScriptFolder+"' doesn't exist! Please run subprogram 'scripts' beforehand.",file=sys.stderr)
        
    jobSubmissionTool    = general["job_submission"]
    for file in os.listdir(exahypeRoot + "/" + outputPath + "/" + scriptsFolder):
        if file.endswith(".job"):
            command=jobSubmissionTool + " " + exahypeRoot + "/" + outputPath + "/" + scriptsFolder + "/" + file
            print(command)
            process = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
            (output, err) = process.communicate()
            process.wait()

if __name__ == "__main__":
    import sys,os
    import configparser
    import subprocess
    import itertools
    import hashlib
    import json
    
    subprograms = ["build","scripts","submit","cleanBuild", "cleanScripts","cleanAll"]
    
    if haveToPrintHelpMessage(sys.argv):
        print("sample usage:./sweep.py myoptions.ini ("+"|".join(subprograms)+")")
        sys.exit()
    
    configFile = parseArgument(sys.argv,1)
    subprogram = parseArgument(sys.argv,2)
    
    config = configparser.ConfigParser()
    config.optionxform=str
    config.read(configFile)
    
    general = config["general"]
    
    environmentspace = parseEnvironment(config)
    parameterspace   = parseParameters(config)
    
    # select subprogram
    if subprogram == "clean":
        clean(general)
    elif subprogram == "cleanBuild":
        clean(general,"build")
    elif subprogram == "cleanScripts":
        clean(general,"scripts")
    elif subprogram == "build":
        build(general,environmentspace,parameterspace)
    elif subprogram == "scripts":
        generateScripts(general,environmentspace,parameterspace)
    elif subprogram == "submit":
        submitJobs(general)
