#!/user/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import re
import os
import csv

import operator

'''
.. module:: extractunigridruntimes
  :platform: Unix, Windows, Mac
  :synopsis: Extract accumulated adapter runtimes from CSV files created
             by the extractadaptertimes.py script.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Extract accumulated adapter runtimes from CSV files created
           by the extractadaptertimes.py script.
'''

def extract_runtimes(table,legacy):
    '''
    Args:
      table (str):
         path of a CSV table created with extractadaptertimes.py
      legacy (bool)
         If true, use old adapter names: ADERDGTimeStep,... Not the new ones: FusedTimeStep,...
    '''
    nonfusedAdapters = [ 
                          "NeighbourDataMerging",
                          "SolutionUpdate",
                          "Prediction" # substract one; exclude initialisation
                       ]
    fusedAdapters    = [ 
                        "FusedTimeStep",
                        "Prediction"
                       ]
    if legacy:
        nonfusedAdapters = [ 
                             "NeighbourDataMerging",
                             "SolutionUpdateAndTimeStepSizeComputation",
                             "Prediction"  # substract one; exclude initialisation
                           ]
        fusedAdapters    = [ 
                            "ADERDGTimeStep",
                            # "PredictionAndFusedTimeSteppingInitialisation", # exclude initialisation
                            "Prediction"
                           ]
    
    datafile    = open(table, 'r')
    header      = next(datafile).strip()
    reader      = csv.reader(datafile,delimiter='&')
    
    # assume: row = [identifier,order,cc,kernels,algorithm,adapter,nodes,tasks,cores,mode,iterations,total user time, total cpu time]
    sorted_data = sorted(reader, key=lambda x: (x[0],x[1],x[2],x[3],x[4],int(x[6]),int(x[7]),int(x[8]),x[5]))
    
    def consider_adapter(row):
        return (row[5].strip() in nonfusedAdapters or 
                row[5].strip() in fusedAdapters)
               
    def only_adapter_changed(row,last_row):
        return (row[0] == last_row[0] and 
                row[1] == last_row[1] and 
                row[2] == last_row[2] and 
                row[3] == last_row[3] and 
                row[4] == last_row[4] and
                # row[5] == last_row[5] and # adapter
                row[6] == last_row[6] and
                row[7] == last_row[7] and
                row[8] == last_row[8])
    
    # init ( we do the first row twice; a little inconsistent )
    last_row = sorted_data[0]
    result = last_row[:5] + last_row[6:]
    result[9]  = '0 '  # iterations
    result[10] = '0 '  # normalised user time
    result[11] = '0'   # normalised fused time step time
    result.append(0)   # predictor reruns
    result.append('0.')# time per predictor rerun
    
    # loop
    header = ["Mesh","Order","CC","Kernels","Algorithm","Nodes","Tasks (per Node)","Cores (per Task)","Mode","Iterations","Total Runtime (Per Iteration)","FusedTimeStep (Per Iteration)","Reruns", "Predictor Rerun (Per Iteration)"]
    filename = table.replace('.csv','.runtimes.csv')
    with open(filename, 'w') as datafile:
        writer = csv.writer(datafile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header)
      
        for row in sorted_data:
            if only_adapter_changed(row,last_row):
                if row[4]=='fused':
                    if row[5]=='FusedTimeStep' or row[5]=='ADERDGTimeStep': # comes before Prediction...
                        result[9]  = row[10] # iterations
                        timesteps = result[9]
                        result[11] = str( float(row[11]) / float(timesteps) ) # normalised user time
                    else:
                        if row[5]=='Prediction':
                            if legacy is False:
                                iterations = float(row[10]) # one more iteration than reruns
                                result[12] = iterations-1
                                result[13] = str ( float(row[11]) / float(iterations) )
                            
                            timesteps = result[9]
                            result[10] = str( float(result[11]) + float(result[12]) * float(result[13]) / float(timesteps) ) # normalised user time
                if row[4]=='nonfused':
                    if row[5]=='NeighbourDataMerging': # comes before Prediction,SolutionUpdate
                        result[9]  = row[10] # iterations
                    if row[5] in nonfusedAdapters:
                        iterations = row[10]
                        result[10] = str( float(result[10]) + float(row[11]) / float(iterations) ) # normalised user time
            else:
                writer.writerow(result)
              
                result = row[:5] + row[6:]
                iterations = row[10]
                result[9]  = '-1'
                result[10] = '0 '  # normalised total runtime
                result[11] = '0'   # normalised fused time step time
                result.append('0') # predictor reruns
                result.append('0.')# time per predictor rerun
            
            last_row = row
        
        # finalise (have to write a last time)
        writer.writerow(result)
    
    # assumes sorted table
#    for algorithm in [ 'fused', 'nonfused' ]:
#        for adapter in adapters[algorithm]:
#            print(algorithm)
#            print(adapter)
#            
#            try:
#                
#                  print(row)
#            except IndexError:
#                print("Line could not be parsed. Skip.")
#            except:
#                print("Something else")
    
    datafile.close() 


#def sort_table(filename):
#    '''
#    Sorts the rows of the file according to nodes,tasks,cores,
#    See: https://stackoverflow.com/a/17109098
#    '''
#    datafile    = open(filename, 'r')
#    header      = next(datafile).strip()
#    reader      = csv.reader(datafile,delimiter='&')
#    # row = [order,cc,kernels,algorithm,nodes,tasks,cores,mode]
#    sorted_data = sorted(reader, key=lambda x: (x[0],int(x[1]),x[2],x[3],x[4],int(x[5]),int(x[6]),int(x[7])))
#    datafile.close() 
# 
#    with open(filename, 'w') as datafile:
#        writer = csv.writer(datafile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
#        writer.writerow(header.split('&'))
#        writer.writerows(sorted_data)

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'extrat_table' above.
help = '''
Extract accumulated adapter runtimes from CSV tables created 
by the extractadaptertimes.py script.
Write the result to another CSV file with name
<table>.runtimes.csv

NOTE: 
This script is solely useful for computing
the runtimes of the fused and nonfused ADER-DG
algorithms on uniform meshes.
Does ignore the mesh setup times.

NOTE: Assumes no plotting was performed during
the simulations!

NOTE: Does not consider AMR or Limiting.
Parsed solve times could be completely wrong!

NOTE: Does not sanitise the input. Be careful that
there are no whitelines in the table file!

\n\n
Sample usage:\n
python extractruntimes.py -table Euler_ADERDG-no-output.csv --legacy
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-table',required=True,help="Directory containing the Peano output files.")
parser.add_argument('--legacy',dest='legacy',required=False,action='store_true',help="Use old adapter names: ADERDGTimeStep,... Not the new ones: FusedTimeStep,...")
parser.set_defaults(legacy=False)

args     = parser.parse_args();

table    = args.table
legacy   = args.legacy

extract_runtimes(table,legacy)
print("created table:")
print(table.replace('.csv','.runtimes.csv'))
