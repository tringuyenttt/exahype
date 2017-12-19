#!/user/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import sys
import re
import os
import csv

import operator

def column(matrix, i):
    return [row[i] for row in matrix]
    

def writeTable(data,header,filename):
    with open(filename, 'w') as datafile:
        writer = csv.writer(datafile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header)
        
        for row in data:
            writer.writerow(row)
    
def extractSingleCoreRuntimes(table):
    '''
    assume: inputRow = [order,adapter,arch,optimisation,run,iteratins,user time, runtime]
    example:
       0           1             2        3      4     5         6        7
    ['6', 'MeshRefinement', 'noarch', 'O3-vec', '1', '12', '15.2857', '15.25'] 
    
    Args:
      table (str):
         path of a CSV table created with extractadaptertimes.py
    '''
    nonfusedAdapters = [ 
                         "BroadcastGlobalDataAndMergeNeighbourMessages",
                         "SolutionUpdate",
                         "Prediction" 
                       ]
    
    datafile = open(table, 'r')
    next(datafile) # skip header
    
    data = list(csv.reader(datafile,delimiter='&'))
    datafile.close()
    
    orders        = set(column(data,0))
    architectures = set(column(data,2))
    optimisations = set(column(data,3))
    
    # processing data
    result = []
    for order in orders:
        for architecture in architectures:
            for optimisation in optimisations:
                averageRuntimePerIteration = 0.0
                minIterations              = sys.maxsize
                for adapter in nonfusedAdapters:
                    filtered = list(filter(lambda x: x[0]==order and x[1]==adapter and x[2]==architecture and x[3]==optimisation, data))
                    runs     = len(filtered)
                    
                    iterations = 0
                    if runs > 0:
                        iterations = int ( filtered[0][5] )
                        
                        averageAdapterTime = 0.0
                        for r in range(0,runs):
                            averageAdapterTime += float ( filtered[r][6] )
                        averageAdapterTime /= runs
                        
                        averageRuntimePerIteration += averageAdapterTime / iterations
                    minIterations = min(minIterations,iterations)
                    
                row = [order,architecture,optimisation,minIterations,averageRuntimePerIteration]
                result.append(row)
    
    header = ["Order","Architecture","Optimisation","Iterations","Total Runtime (Per Iteration)"]
    filename = table.replace('.csv','.runtimes.csv')
    writeTable(result,header,filename)
'''
.. module:: extractunigridruntimes
  :platform: Unix, Windows, Mac
  :synopsis: Extract accumulated adapter runtimes from CSV files created
             by the extractadaptertimes.py script.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Extract accumulated adapter runtimes from CSV files created
           by the extractadaptertimes.py script.
'''

def extractRuntimes(table):
    '''
    Args:
      table (str):
         path of a CSV table created with extractadaptertimes.py
      single-core (bool)
         If true, use old adapter names: ADERDGTimeStep,... Not the new ones: FusedTimeStep,...
    '''
    nonfusedAdapters = [ 
                          "BroadcastGlobalDataAndMergeNeighbourMessages",
                          "SolutionUpdate",
                          "Prediction" # substract one; exclude initialisation
                       ]
    fusedAdapters    = [ 
                        "FusedTimeStep",
                        "PredictionRerun"
                       ]
    
    datafile    = open(table, 'r')
    header      = next(datafile).strip()
    reader      = csv.reader(datafile,delimiter='&')
    
    # assume: row = [identifier,order,cc,kernels,algorithm,adapter,nodes,tasks,cores,mode,iterations,total user time, total cpu time]
    sorted_data = sorted(reader, key=lambda x: (x[0],x[1],x[2],x[3],x[4],int(x[6]),int(x[7]),int(x[8]),x[5]))
               
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
                            if single-core is False:
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
    
    datafile.close()

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
python extractruntimes.py -table Elastic3D-no-output.csv --single-core
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-table',required=True,help="Directory containing the Peano output files.")
parser.add_argument('--single-core',dest='singlecore',required=False,action='store_true',help="Read a single-core table which is structured differently than the other tables.")
parser.set_defaults(singlecore=False)

args = parser.parse_args();

table      = args.table
singlecore = args.singlecore

if singlecore is False:
  extractRuntimes(table)
else:
  extractSingleCoreRuntimes(table)  

print("created table:")
print(table.replace('.csv','.runtimes.csv'))
