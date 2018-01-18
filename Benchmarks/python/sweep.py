#!/usr/bin/python3
'''
.. module:: sweep
  :platform: Unix, Windows, Mac
  :synopsis: Generate benchmark suites for ExaHyPE.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>, 

:synopsis: Generate benchmark suites for ExaHyPE.
'''
class Subprogram:
    def run(self):
        print('Hello, how are you?')

def extractSubprogram(argv):
    if len(argv)>1:
        return argv[1]
    else:
        return None

def extractOptionsFile(argv):
    if len(argv)>2:
        return argv[2]
    else:
        return None

def haveToPrintHelpMessage(argv):
    result = extractSubprogram(argv) ==None or \
             extractOptionsFile(argv)==None
    print(result)
    for arg in argv:
        result = result or ( arg=='-help' or arg=='-h' )
    print(result)
    return result

if __name__ == '__main__':
    import sys,os
    
    # for debugging
    #print("Number of arguments: "+str(len(sys.argv)))
    #print("The arguments are: "+str(sys.argv))
    #print("The subprogram is: "+str(sys.argv[1]))
    #print("The options file is: "+extractOptionsFile(sys.argv))
    os.environ["MY_TEST_ENV"] = "1" # child processes (build processes) will inherit this environment variable
    print(os.environ)
    
    if haveToPrintHelpMessage(sys.argv):
        print("sample usage: python sweep (setup|build|generate) options.sweep")
        print("sample usage: python sweep (setup|build|generate) options.sweep -order 3 5 7")
    
