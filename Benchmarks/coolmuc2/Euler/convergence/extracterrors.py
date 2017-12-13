#!/user/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import re
import os
import csv

import operator

'''
.. module:: extracterrors
  :platform: Unix, Windows, Mac
  :synopsis: Extracts errors from Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Extracts errors from Peano output files with specific file naming pattern.
'''


def parse_errors(file_path):
    """
    Reads a single Peano output file and parses errors
    
    Args:
       file_path (str):
          Path to the Peano output file.

    Returns:
       A dict holding for each of the variables a dict with a list for absErrorL1,~L2,~LInf,relErrorL1,~L2,~LInf
    """
    result = { }
    errors = ["absErrorL1","absErrorL2","absErrorLInf","relErrorL1","relErrorL2","relErrorLInf"]

    try:
        file_handle=open(file_path)
        
        for line in file_handle:
            if "coarsest mesh size was chosen as " in line:
                match = re.search('chosen as ([0-9]+\.[0-9]+)',line)
                hMax = float(match.group(1).strip())
                result["h"] = hMax
            i = 0
            if "t_eval" in line:
                result["t"][i] = float(line.split(':')[-1].strip())
            variables = 0;
            if i==0 and "variable" in line
                variables = len(split(':')[1].split(','))
                i+=1             
            for variable in range(0,variables):
                result["var"+variable] = {}
                for error in errors:
                    result["var"+variable][error] = {}

            for error in errors:
                if error in line:
                    segments = split(':')[1].split(',')
                    for variable in range(0,variables):
                        result["var"+variable][error][i] = float(segments[variable]);
    except:
        print ("Error: Could not process file '%s'!\n" % (file_path))
        raise
    return result

def extract_table(root_dir,prefix):
    '''
    Extracts performance metrics from Peano output files with specific file naming pattern.
   
    Args:
      root_dir (str):
         Directory containing the Peano output files.
      prefix (str):
         Prefix of the files - usually the date of the test and an identifier for the test.
    '''
    errors = ["absErrorL1","absErrorL2","absErrorLInf","relErrorL1","relErrorL2","relErrorLInf"]
    
    # collect filenames
    with open(prefix+'.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for filename in os.listdir(root_dir):
            if filename.endswith(".out") and filename.startswith(prefix):
                # sample: EulerADERDG-gen-regular-0-fused-p3.out
                match = re.search('^.+-.+-([a-z]+)-.+-([A-Za-z]+)-p([0-9]+)',filename)
                opt   = match.group(1)
                fused = match.group(2)
                order = match.group(3)
                    
                result = parse_errors(filename) 
                
                h    = result["h"]
                N    = len(result["t"])
                Nvar = len(result)-2
                for variable in range(0,Nvar):
                    for i in range(0,N):
                        csvwriter.writerow([
                            result["h"], variable, result["t"][i],
                            result["var"+variable]["absErrorL1"]  [i],
                            result["var"+variable]["absErrorL2"]  [i],
                            result["var"+variable]["absErrorLInf"][i],
                            result["var"+variable]["relErrorL1"]  [i],
                            result["var"+variable]["relErrorL2"]  [i],
                            result["var"+variable]["relErrorLInf"][i] 
                        ]);

def sort_table(filename):
    '''
    Sorts the rows of the file according to nodes,tasks,cores,adapter name.
    See: https://stackoverflow.com/a/17109098
    '''
    datafile    = open(filename, 'r')
    reader      = csv.reader(datafile,delimiter='&')
    sorted_data = sorted(reader, key=lambda x: (x[3],int(x[0]),int(x[1]),int(x[2])))
    datafile.close() 
 
    with open(filename, 'w') as datafile:
        writer = csv.writer(datafile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(sorted_data)

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'extrat_table' above.
help = '''
Extract performance metrics from Peano output files with specific file naming pattern
and write them to a csv file with name
<prefix>.csv

\n\n
Sample usage:\n
python3 extracterrors.py -path \'results/' -prefix \'Euler-ADERDG-regular-0-fused-p3\'
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-path',required=True,help="Directory containing the Peano output files.")
parser.add_argument('-prefix',required=True,help="Prefix of the Peano output files.")

args     = parser.parse_args();

root_dir = args.path
prefix   = args.prefix

extract_table(root_dir,prefix)
