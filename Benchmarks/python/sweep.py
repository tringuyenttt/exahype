#!/usr/bin/env python3
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
    environmentSpace = {}
    if "environment" in config and len(config["environment"].keys()):
        for key, value in config["environment"].items():
            environmentSpace[key] = [x.strip() for x in value.split(",")]
        if "SHAREDMEM" not in environmentSpace:
            print("ERROR: 'SHAREDMEM' missing in section 'environment'.",file=sys.stderr)
            sys.exit()
    else:
        print("ERROR: Section 'environment' is missing or empty! (Must contain at least 'SHAREDMEM'.)",file=sys.stderr)
        sys.exit()
    
    return environmentSpace


def parseParameters(config):
    """
    Parse the parameters section.
    """
    parameterSpace = {}
    if "parameters" in config and len(config["parameters"].keys()):
        for key, value in config["parameters"].items():
            parameterSpace[key] = [x.strip() for x in value.split(",")]
            
        if "order" not in parameterSpace:
            print("ERROR: 'order' missing in section 'parameters'.",file=sys.stderr)
            sys.exit()
        elif "dimension" not in parameterSpace:
            print("ERROR: 'dimension' missing in section 'parameters'.",file=sys.stderr)
            sys.exit()
    else:
        print("ERROR: Section 'parameters' is missing or empty! (Must contain at least 'dimension' and 'order'.)",file=sys.stderr)
        sys.exit()
    
    return parameterSpace

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
    
    return hashlib.md5(chain.encode()).hexdigest()

usedEnvironmentDictHashes = []
usedParameterDictHashes   = []

def hashEnvironmentDict(dictionary):
    """
    Hash a dictionary.
    """
    result = hashDictionary(dictionary)
    
    if result not in usedParameterDictHashes:
        usedParameterDictHashes.append(result)
    else:
        print("ERROR: Hash conflict for environment option sets. Try to resolve it by adding another (dummy) environment variable.",file=sys.stderr)
        sys.exit()
    return result
    
def hashParameterDict(dictionary):
    """
    Hash a dictionary.
    """
    result = hashDictionary(dictionary)
    
    if result not in usedParameterDictHashes:
        usedParameterDictHashes.append(result)
    else:
        print("ERROR: Hash conflict for parameter option sets. Try to resolve it by adding another (dummy) parameter.",file=sys.stderr)
        sys.exit()
    return result

def clean(subFolder=""):
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

def build(buildOnlyMissing=False):
    """
    Build the executables.
    """
    print("currently loaded modules:")
    subprocess.call("module list",shell=True)
    print("")
    print("ExaHyPE build environment:")
    exahypeEnv = ["COMPILER", "MODE", "SHAREDMEM", "DISTRIBUTEDMEM", "EXAHYPE_CC", "PROJECT_C_FLAGS", "PROJECT_L_FLAGS", "COMPILER_CFLAGS", "COMPILER_LFLAGS", "FCOMPILER_CFLAGS", "FCOMPILER_LFLAGS"]
    for variable in exahypeEnv:
        if variable in os.environ:
            print(variable+"="+os.environ[variable])
    print("")
    
    templateFileName = general["spec_template"]
    exahypeRoot      = general["exahype_root"]
    outputPath       = general["output_path"]
    projectPath      = general["project_path"]
    projectName      = general["project_name"]
    
    templateBody = None
    with open(exahypeRoot+"/"+templateFileName, "r") as templateFile:
        templateBody=templateFile.read()
    
    if templateBody!=None:
        buildFolderPath = exahypeRoot+"/"+outputPath+"/"+buildFolder
        
        if not os.path.exists(buildFolderPath):
            print("create directory:"+buildFolderPath)
            os.makedirs(buildFolderPath)
        
        dimensions = parameterSpace["dimension"]
        orders     = parameterSpace["order"]
        buildParameterDict = list(dictProduct(parameterSpace))[0]
        
        firstIteration = True
        newExecutables=0
        for environmentDict in dictProduct(environmentSpace):
            for key,value in environmentDict.items():
                os.environ[key]=value
            
            for dimension in dimensions:
                if not firstIteration:
                    command = "make clean"
                    print(command)
                    process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    (output, err) = process.communicate()
                    process.wait()
                
                for order in orders:
                    oldExecutable = exahypeRoot + "/" + projectPath+"/ExaHyPE-"+projectName
                    newExecutable = buildFolderPath + "/ExaHyPE-"+projectName+"-"+hashEnvironmentDict(environmentDict)+"-d" + dimension + "-p" + order
                    
                    if not os.path.exists(newExecutable) or not buildOnlyMissing:
                        buildParameterDict["dimension"]=dimension
                        buildParameterDict["order"]    =order
                        
                        buildSpecFileBody = renderSpecificationFile(templateBody,buildParameterDict,"1","1")
                        
                        buildspecFilePath = outputPath+"/"+buildFolder+"/"+projectName+"-d"+dimension+"-p"+order+".exahype"
                        
                        with open(exahypeRoot + "/" + buildspecFilePath, "w") as buildSpecificationFile:
                            buildSpecificationFile.write(buildSpecFileBody)
                        
                        print("building executable for " + \
                          "environment="+str(environmentDict) + \
                          ", dimension="+dimension + \
                          ", order="+order,file=sys.stderr)
                        # run toolkit
                        toolkitCommand = "(cd "+exahypeRoot+" && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive "+buildspecFilePath+")"
                        print(toolkitCommand,end="",flush=True)
                        process = subprocess.Popen([toolkitCommand], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                        (output, toolkitErr) = process.communicate()
                        process.wait()
                        if "setup build environment ... ok" in str(output):
                            print(" [OK]")
                        else:
                            print(" [FAILED]")
                            print("toolkit errors/warnings=\n"+toolkitErr.decode('UTF-8'),file=sys.stderr)
                            sys.exit()
                        
                        if firstIteration:
                            command = "make clean"
                            print(command)
                            process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                            (output, err) = process.communicate()
                            process.wait()
                            firstIteration = False
                        else: # clean application folder only
                            command = "rm -r *.o cipofiles.mk cfiles.mk ffiles.mk kernels"
                            print(command)
                            process = subprocess.Popen(["make clean"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                            process.communicate()
                            process.wait()
                        
                        # call make
                        make_threads=general["make_threads"]
                        makeCommand="make -j"+make_threads
                        print(makeCommand,end="",flush=True)
                        process = subprocess.Popen([makeCommand], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                        (output, makeErr) = process.communicate()
                        process.wait()
                        if "build of ExaHyPE successful" in str(output):
                            print(" [OK]")
                        else:
                            print(" [FAILED]")
                            print("make errors/warnings=\n"+makeErr.decode('UTF-8'),file=sys.stderr)
                            sys.exit()
                            
                        moveCommand   = "mv "+oldExecutable+" "+newExecutable
                        print(moveCommand)
                        subprocess.call(moveCommand,shell=True)
                        print("--------------------------------------------------------------------------------")
                        print("toolkit errors/warnings=\n"+toolkitErr.decode('UTF-8'),file=sys.stderr)
                        print("make errors/warnings=\n"+makeErr.decode('UTF-8'),file=sys.stderr)
                        print("--------------------------------------------------------------------------------")
                        newExectutables+=1
                    else:
                        print("skipped building of '"+newExecutables+"' as it already exists.")
        print("built "+newExecutables+" executables")
    else:
        print("ERROR: couldn\'t open template file: "+templateFileName,file=sys.stderr)

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
                    jobName,jobFilePath,outputFileName,errorFileName,appName,specFilePath,
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
    
    context["job_file"]    = jobFilePath
    context["app"]        = appName
    context["spec_file"]  = specFilePath
    
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

def verifyAllExecutablesExist(justWarn=False):
    exahypeRoot = general["exahype_root"]
    outputPath  = general["output_path"]
    projectName = general["project_name"]
    
    dimensions = parameterSpace["dimension"]
    orders     = parameterSpace["order"]
    
    messageType = "ERROR"
    if justWarn:
      messageType = "WARNING"
    
    buildFolderPath = exahypeRoot+"/"+outputPath+"/"+buildFolder
    if not justWarn and not os.path.exists(buildFolderPath):
        print("ERROR: build folder '"+buildFolderPath+"' doesn't exist! Please run subprogram 'build' beforehand.",file=sys.stderr)
        sys.exit()
    
    allExecutablesExist = True
    for environmentDict in dictProduct(environmentSpace):
        environmentDictHash = hashEnvironmentDict(environmentDict)
        for dimension in dimensions:
            for order in orders:
                appName   = exahypeRoot + "/" + outputPath+"/"+buildFolder + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-d" + dimension + "-p" + order
                
                if not os.path.exists(appName):
                    allExecutablesExist = False
                    print(messageType+ ": application for " + \
                          "environment="+str(environmentDict) + \
                          ", dimension="+dimension + \
                          ", order="+order + \
                          " does not exist! ('"+appName+"')",file=sys.stderr)
    
    if not justWarn and not allExecutablesExist:
        print("ERROR: subprogram failed as not all executables exist. Please adopt your options file according to the error messages.\n" + \
              "       Then rerun the 'build' subprogram.",file=sys.stderr)
        sys.exit()

def generateScripts():
    """
    Generate spec files and job scripts.
    """
    exahypeRoot   = general["exahype_root"]
    outputPath    = general["output_path"]
    projectName   = general["project_name"]
    
    jobs       = config["jobs"]
    nodeCounts = [x.strip() for x in jobs["nodes"].split(",")]
    taskCounts = [x.strip() for x in jobs["tasks"].split(",")]
    coreCounts = parseCores(jobs);
    runs       = int(jobs["runs"])
    
    specFileTemplatePath  = exahypeRoot+"/"+general["spec_template"]
    jobScriptTemplatePath = exahypeRoot+"/"+general["job_template"]
    
    specFileTemplate  = None
    with open(specFileTemplatePath, "r") as templateFile:
        specFileTemplate=templateFile.read()
    if specFileTemplate is None:
        print("ERROR: couldn\'t open template file: "+specFileTemplatePath,file=sys.stderr)
        sys.exit()
        
    jobScriptTemplate = None
    with open(jobScriptTemplatePath, "r") as templateFile:
        jobScriptTemplate=templateFile.read()
    if jobScriptTemplate is None:
        print("ERROR: couldn\'t open template file: "+jobScriptTemplatePath,file=sys.stderr)
        sys.exit()
        
    if not os.path.exists(exahypeRoot+"/"+outputPath+"/"+scriptsFolder):
        os.makedirs(exahypeRoot+"/"+outputPath+"/"+scriptsFolder)
    if not os.path.exists(exahypeRoot+"/"+outputPath+"/"+resultsFolder):
        os.makedirs(exahypeRoot+"/"+outputPath+"/"+resultsFolder)
    
    # spec files
    specFile=0
    for parameterDict in dictProduct(parameterSpace):
        parameterDictHash = hashParameterDict(parameterDict)
        
        for tasks in taskCounts:
            for cores in coreCounts:
              specFileBody = renderSpecificationFile(specFileTemplate,parameterDict,tasks,cores)
              
              specFilePath = exahypeRoot + "/" + outputPath + "/" + scriptsFolder + "/" + projectName + "-" + parameterDictHash + "-t"+tasks+"-c"+cores+".exahype"
              
              with open(specFilePath, "w") as specFile:
                  specFile.write(specFileBody)
              specFiles+=1
    
    print("generated "+specFiles+" specification files.")
    
    # check if required executables exist
    verifyAllExecutablesExist(True)
    
    # generate job scrips
    jobScripts = 0
    for run in range(0,runs):
        for nodes in nodeCounts:
            for tasks in taskCounts:
                for cores in coreCounts:
                    for environmentDict in dictProduct(environmentSpace):
                        environmentDictHash = hashDictionary(environmentDict)
                        
                        for parameterDict in dictProduct(parameterSpace):
                            parameterDictHash = hashDictionary(parameterDict)
                            
                            dimension = parameterDict["dimension"]
                            order     = parameterDict["order"]
                            appName   = exahypeRoot + "/" + outputPath+"/"+buildFolder + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-d" + dimension + "-p" + order
                            
                            
                            specFilePath = exahypeRoot + "/" + outputPath + "/" + scriptsFolder + "/" + projectName + "-" + \
                                           parameterDictHash + "-t"+tasks+"-c"+cores+".exahype"
                            
                            jobName        = projectName + "-" + environmentDictHash + "-" + parameterDictHash + \
                                             "-n" + nodes + "-t"+tasks+"-c"+cores+"-r"+str(run)
                            jobFilePrefix  = exahypeRoot + "/" + outputPath + "/" + scriptsFolder + "/" + jobName
                            jobFilePath    = jobFilePrefix + ".job"
                            outputFileName = exahypeRoot + "/" + outputPath + "/" + resultsFolder + "/" + jobName + ".out"
                            errorFileName  = exahypeRoot + "/" + outputPath + "/" + resultsFolder + "/" + jobName + ".err"
                            
                            jobScriptBody = renderJobScript(jobScriptTemplate,environmentDict,parameterDict,jobs,
                                                            jobName,jobFilePath,outputFileName,errorFileName,appName,specFilePath,
                                                            nodes,tasks,cores,run)
                            with open(jobFilePath, "w") as jobFile:
                                jobFile.write(jobScriptBody)
                            
                            jobScripts+=1

    print("generated "+jobScripts+" job scripts")

                             
def verifyAllJobScriptsExist():
    """
    Verify that all job scripts exist.
    """
    exahypeRoot          = general["exahype_root"]
    outputPath           = general["output_path"]
    projectName          = general["project_name"]
    jobSubmissionTool    = general["job_submission"]
    
    jobs       = config["jobs"]
    nodeCounts = [x.strip() for x in jobs["nodes"].split(",")]
    taskCounts = [x.strip() for x in jobs["tasks"].split(",")]
    coreCounts = parseCores(jobs);
    runs       = int(jobs["runs"])
    
    scriptFolderPath = exahypeRoot+"/"+outputPath+"/"+scriptsFolder
    if not os.path.exists(scriptFolderPath):
        print("ERROR: job script folder '"+scriptFolderPath+"' doesn't exist! Please run subprogram 'scripts' beforehand.",file=sys.stderr)
        sys.exit()
    
    allJobScriptsExist = True
    for run in range(0,runs):
        for nodes in nodeCounts:
            for tasks in taskCounts:
                for cores in coreCounts:
                    for environmentDict in dictProduct(environmentSpace):
                        environmentDictHash = hashDictionary(environmentDict)
                        
                        for parameterDict in dictProduct(parameterSpace):
                            parameterDictHash = hashDictionary(parameterDict)
                            
                            dimension = parameterDict["dimension"]
                            order     = parameterDict["order"]
                            
                            jobName        = projectName + "-" + environmentDictHash + "-" + parameterDictHash + \
                                             "-n" + nodes + "-t"+tasks+"-c"+cores+"-r"+str(run)
                            jobFilePrefix  = exahypeRoot + "/" + outputPath + "/" + scriptsFolder + "/" + jobName
                            jobFilePath    = jobFilePrefix + ".job"
                            
                            if not os.path.exists(jobFilePath):
                                allJobScriptsExist = False
                                print("ERROR: job script for " + \
                                      "environment="+str(environmentDict)+ \
                                      ", parameters="+str(parameterDict) + \
                                      ", nodes="+nodes + \
                                      ", tasks="+tasks + \
                                      ", cores="+cores + \
                                      ", run="+str(run) + \
                                      " does not exist! ('"+jobFilePath+"')",file=sys.stderr)
    if not allJobScriptsExist:
        print("ERROR: subprogram failed! Please adopt your sweep options file according to the error messages.\n" + \
              "       Then rerun the 'scripts' subprogram.")
        sys.exit()

def verifyAllSpecFilesExist():
    """
    Verify that all ExaHyPE specification files exist.
    """
    exahypeRoot          = general["exahype_root"]
    outputPath           = general["output_path"]
    projectName          = general["project_name"]
    jobSubmissionTool    = general["job_submission"]
    
    jobs       = config["jobs"]
    taskCounts = [x.strip() for x in jobs["tasks"].split(",")]
    coreCounts = parseCores(jobs);
    
    scriptFolderPath = exahypeRoot+"/"+outputPath+"/"+scriptsFolder
    if not os.path.exists(scriptFolderPath):
        print("ERROR: job script folder '"+scriptFolderPath+"' doesn't exist! Please run subprogram 'scripts' beforehand.",file=sys.stderr)
        sys.exit()
    
    allSpecFilesExist = True
    for parameterDict in dictProduct(parameterSpace):
        parameterDictHash = hashParameterDict(parameterDict)
        
        for tasks in taskCounts:
            for cores in coreCounts:
              specFilePath = exahypeRoot + "/" + outputPath + "/" + scriptsFolder + "/" + projectName + "-" + parameterDictHash + "-t"+tasks+"-c"+cores+".exahype"
              
              if not os.path.exists(specFilePath):
                  allSpecFilesExist = False
                  print("ERROR: specification file for \n" + \
                        "parameters="+str(parameterDict) + \
                        ", tasks="+tasks + \
                        ", cores="+cores + \
                        " does not exist! ('"+specFilePath+"')",file=sys.stderr)
    
    if not allSpecFilesExist:
        print("ERROR: subprogram failed! Please adopt your sweep options file according to the error messages.\n" + \
              "       Then rerun the 'scripts' subprogram.")
        sys.exit()

def hashSweep(jobs,enviromentSpace,parameterSpace):
    nodeCounts = [x.strip() for x in jobs["nodes"].split(",")]
    taskCounts = [x.strip() for x in jobs["tasks"].split(",")]
    coreCounts = parseCores(jobs);
    runs       = jobs["runs"]
    
    chain = ""
    for value in nodeCounts:
        chain += value+";"
    for value in taskCounts:
        chain += value+";"
    for value in coreCounts:
        chain += value+";"
    chain += runs+";"
    
    for environmentDict in dictProduct(environmentSpace):
        chain += hashDictionary(environmentDict)
    for parameterDict in dictProduct(parameterSpace):
        chain += hashDictionary(parameterDict)
        
    return hashlib.md5(chain.encode()).hexdigest()

def extractJobId(processOutput):
    jobId = "unknown"
    if "Submitted batch job " in processOutput:
        jobId = processOutput.split(" ")[-1]
    return jobId

def submitJobs():
    exahypeRoot          = general["exahype_root"]
    outputPath           = general["output_path"]
    projectName          = general["project_name"]
    jobSubmissionTool    = general["job_submission"]
    
    jobs       = config["jobs"]
    nodeCounts = [x.strip() for x in jobs["nodes"].split(",")]
    taskCounts = [x.strip() for x in jobs["tasks"].split(",")]
    coreCounts = parseCores(jobs);
    runs       = int(jobs["runs"])
    
    # verify everything is fine
    verifyAllExecutablesExist()
    verifyAllJobScriptsExist()
    verifyAllSpecFilesExist()
    
    # loop over job scrips
    jobIds = []
    for run in range(0,runs):
        for nodes in nodeCounts:
            for tasks in taskCounts:
                for cores in coreCounts:
                    for environmentDict in dictProduct(environmentSpace):
                        environmentDictHash = hashDictionary(environmentDict)
                        
                        for parameterDict in dictProduct(parameterSpace):
                            parameterDictHash = hashDictionary(parameterDict)
                            
                            dimension = parameterDict["dimension"]
                            order     = parameterDict["order"]
                            
                            jobName        = projectName + "-" + environmentDictHash + "-" + parameterDictHash + \
                                             "-n" + nodes + "-t"+tasks+"-c"+cores+"-r"+str(run)
                            jobFilePrefix  = exahypeRoot + "/" + outputPath + "/" + scriptsFolder + "/" + jobName
                            jobFilePath    = jobFilePrefix + ".job"
                            
                            command=jobSubmissionTool + " " + jobFilePath
                            print(command)
                            process = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
                            (output, err) = process.communicate()
                            process.wait()
                            jobIds.append(extractJobId(output.decode("UTF_8"))
                            
    submittedJobsPath = exahypeRoot + "/" + outputPath + "/" + \
                        hashSweep(jobs,environmentSpace,parameterSpace) + ".submitted"
    
    print(jobIds)
    with open(submittedJobsPath, "w") as submittedJobsFile:
        submittedJobsFile.write(json.dumps(jobIds))
    
    print("submitted "+str(len(jobIds))+" jobs")
    print("job ids are memorised in: "+submittedJobsPath)

def cancelJobs():
    exahypeRoot         = general["exahype_root"]
    outputPath          = general["output_path"]
    jobCancellationTool = general["job_cancellation"]
 
    submittedJobsPath = exahypeRoot + "/" + outputPath + "/" + \
                        hashSweep(jobs,environmentSpace,parameterSpace) + ".submitted"    

    jobIds = None
    with open(submittedJobsPath, "r") as submittedJobsFile:
        jobIds = json.loads(submittedJobsFile.read())
    
    if jobIds==None:
        print("ERROR: couldn't find any submitted jobs for current sweep ('"+submittedJobsPath+"').")
        sys.exit()
     
    for jobId in jobIds:
        command = jobCancellationTool + " " + jobId
        print(command)
        subprocess.call(command,shell=True)
    print("cancelled "+str(len(jobIds))+" jobs")    

    command = "rm "+submittedJobsPath
    print(command)
    subprocess.call(command,shell=True)

if __name__ == "__main__":
    import sys,os
    import configparser
    import subprocess
    import itertools
    import hashlib
    import json
    import codecs
    
    subprograms = ["build","buildMissing","scripts","submit","cancel","parseResults","cleanBuild", "cleanScripts","cleanAll"]
    scriptsFolder        = "scripts"
    buildFolder          = "build"
    resultsFolder        = "results"
    
    if haveToPrintHelpMessage(sys.argv):
        print("sample usage:./sweep.py myoptions.ini ("+"|".join(subprograms)+")")
        sys.exit()
    
    configFile = parseArgument(sys.argv,1)
    subprogram = parseArgument(sys.argv,2)
    
    config = configparser.ConfigParser()
    config.optionxform=str
    config.read(configFile)
    
    general          = config["general"]
    jobs             = config["jobs"]
    environmentSpace = parseEnvironment(config)
    parameterSpace   = parseParameters(config)
    
    # select subprogram
    if subprogram == "clean":
        clean(general)
    elif subprogram == "cleanBuild":
        clean("build")
    elif subprogram == "cleanScripts":
        clean("scripts")
    elif subprogram == "build":
        build()
    elif subprogram == "buildMissing":
        build(True)
    elif subprogram == "scripts":
        generateScripts()
    elif subprogram == "submit":
        submitJobs()
    elif subprogram == "cancel":
        cancelJobs()
    elif subprogram == "parseResults":
        print("Not implemented yet!")
