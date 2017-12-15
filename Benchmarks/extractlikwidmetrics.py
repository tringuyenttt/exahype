#!/user/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import re
import os
import csv

import operator

'''
.. module:: writelikwidcsvtable
  :platform: Unix, Windows, Mac
  :synopsis: Contains routines to extract 
             performance metrics from Peano output files with specific file naming pattern.
             The call of the Peano application must be wrapped within likwid-perfctr via: 
             likwid-perfctr -C 0 -g <METRIC> <application>
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Extracts performance metrics from Peano output files with specific file naming pattern.
'''

def extract_likwid_metrics(root_dir,prefix):
    '''
    Extracts performance metrics from Peano output files with specific file naming pattern.
    
    We are currently parsing  "MFLOP/s STAT" (SUM), "Memory bandwidth [MBytes/s]" (SUM),
    "Branch misprediction rate STAT" (Avg),  "L2 miss rate STAT" (Avg).
    
    The output of likwid is different for singlecore and multicore
    runs:
    
    # Singlecore Metrics

    ## FLOPS_DP

    +----------------------+-----------+
    |        Metric        |   Core 1  |
    +----------------------+-----------+
    |  Runtime (RDTSC) [s] |   6.4970  |
    | Runtime unhalted [s] |   7.7641  |
    |      Clock [MHz]     | 2867.1614 |
    |          CPI         |   1.4365  |
    |        MFLOP/s       |  578.7295 |
    |      AVX MFLOP/s     |     0     |
    |    Packed MUOPS/s    |  252.4243 |
    |    Scalar MUOPS/s    |  73.8808  |
    +----------------------+-----------+
    
    # Multicore Metrics
    
    
    ## FLOPS_DP:
    
    +---------------------------+------------+-----------+-----------+-----------+
    |           Metric          |     Sum    |    Min    |    Max    |    Avg    |
    +---------------------------+------------+-----------+-----------+-----------+
    |  Runtime (RDTSC) [s] STAT |   27.7776  |   1.1574  |   1.1574  |   1.1574  |
    | Runtime unhalted [s] STAT |   17.5755  |   0.7025  |   1.1664  |   0.7323  |
    |      Clock [MHz] STAT     | 53287.5604 | 2194.3481 | 2339.6363 | 2220.3150 |
    |          CPI STAT         |   80.7068  |   0.8296  |   3.5650  |   3.3628  |
    |        MFLOP/s STAT       |  3248.7237 |  118.0829 |  532.8151 |  135.3635 |
    |      AVX MFLOP/s STAT     |      0     |     0     |     0     |     0     |
    |    Packed MUOPS/s STAT    |  1416.9944 |  59.0414  |  59.0415  |  59.0414  |
    |    Scalar MUOPS/s STAT    |  414.7347  |     0     |  414.7321 |  17.2806  |
    +---------------------------+------------+-----------+-----------+-----------+
    
    # Counters
    
    Counters have the same singlecore-multicore difference and additionally an extra column:
    
    +-----------------------------------------------+---------+----------------+-------------+---------------+--------------+
    |                     Event                     | Counter |       Sum      |     Min     |      Max      |      Avg     |
    +-----------------------------------------------+---------+----------------+-------------+---------------+--------------+
    |             INSTR_RETIRED_ANY STAT            |  FIXC0  | 11181192786999 | 88372626842 | 1091629701014 | 9.317661e+11 |
    |           CPU_CLK_UNHALTED_CORE STAT          |  FIXC1  |  3312729891579 | 48608952127 |  325502991727 | 2.760608e+11 |
    |           CPU_CLK_UNHALTED_REF STAT           |  FIXC2  |  2930209122522 | 43828197666 |  284406478114 | 2.441841e+11 |
    | FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE STAT |   PMC0  |    389255648   |      0      |    37863481   | 3.243797e+07 |
    |    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  |  623010105225  |    11437    |  66306184117  | 5.191751e+10 | 
    | FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE STAT |   PMC2  |   6888074859   |      0      |   641967988   | 5.740062e+08 |
    +-----------------------------------------------+---------+----------------+-------------+---------------+--------------+

    
    Args:
      root_dir (str):
         Directory containing the Peano output files.
      prefix (str):
         Prefix of the files - usually the date of the test and an identifier for the test.
    '''
    
    metrics    = [
                  "MFLOP/s STAT",
                  "Memory bandwidth [MBytes/s]",
                  "Memory data volume [GBytes] STAT",
                  "Branch misprediction rate STAT",
                  "L2 miss rate STAT" 
                  ]
                  
    columns    = [ "Sum","Min","Max","Avg" ]
    
    # collect filenames
    with open(prefix+'.likwid.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
        # write header
        header = ["Nodes","Tasks (per Node)","Cores (per Task)","CC","Mode"]
        for metric in metrics:
            for column in columns:          
                    header.append(metric+"("+column+")")
        csvwriter.writerow(header)

        # write content
        for filename in os.listdir(root_dir):
            if filename.endswith(".out.likwid") and filename.startswith(prefix):
                print(filename)
                match = re.search(prefix+'-n([0-9]+)-t([0-9]+)-c([0-9]+)-([A-Za-z]+)-([A-Za-z]+)\.out',filename)
                nodes = match.group(1)
                tasks = match.group(2)
                cores = match.group(3)
                mode  = match.group(4)
                cc    = match.group(5)
                    
                measurements = parse_likwid_metrics(root_dir+'/'+filename,metrics,int(cores)==1) 
                
                row = [nodes,tasks,cores,cc,mode]
                for metric in metrics:
                    for column in columns:
                        row.append ( str(measurements[metric][column]) )
                csvwriter.writerow(row)

def parse_likwid_metrics(file_path,metrics,singlecore=False):
    """
    Reads a single Peano output file and parses likwid performance metrics.
    
    Args:
       file_path (str):
          Path to the Peano output file.
       metrics (str[]):
          A list of the metrics we want to read out.
       singlecore (bool):
          Specifies if the run was a singlecore run.

    Returns:
       A dict holding for each of the found metrics a nested dict that holds the following key-value pairs:
          * 'Sum' 
          * 'Avg' 
          * 'Min' 
          * 'Max' 
    """
    columns    = [ "Sum","Min","Max","Avg" ]
    
    result  = { }
    for metric in metrics:
        result[metric] =  { }
        for column in columns:
            result[metric][column] = -1.0

    try:
        file_handle=open(file_path)
        
        for line in file_handle:
            for metric in metrics: 
                if metric in line:
                    segments = line.split('|')
                        
                    if singlecore:
                        #    |     Runtime (RDTSC) [s]    |    6.5219    |
                        value = float(segments[2].strip());
                        result[metric]["Sum"] = value
                        result[metric]["Min"] = value
                        result[metric]["Max"] = value
                        result[metric]["Avg"] = value
                    else:
                        #   |  Runtime (RDTSC) [s] STAT |   27.4632  |   1.1443  |   1.1443  |   1.1443  |
                        result[metric]["Sum"] = float(segments[2].strip());
                        result[metric]["Min"] = float(segments[3].strip());
                        result[metric]["Max"] = float(segments[4].strip());
                        result[metric]["Avg"] = float(segments[5].strip());
    except:
        print ("Error: Could not process file '%s'!\n" % (file_path))
        raise
    return result


def sort_table(filename):
    '''
    Sorts the rows of the file according to nodes,tasks,cores,
    See: https://stackoverflow.com/a/17109098
    '''
    datafile    = open(filename, 'r')
    header      = next(datafile).strip()
    reader      = csv.reader(datafile,delimiter='&')
    sorted_data = sorted(reader, key=lambda x: (int(x[0]),int(x[1]),int(x[2])))
    datafile.close() 
 
    with open(filename, 'w') as datafile:
        writer = csv.writer(datafile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header.split('&'))
        writer.writerows(sorted_data)

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'extrat_table' above.
help = '''
Extract performance metrics from Peano output files with specific file naming pattern
and write them to a csv file with name
<prefix.csv

\n\n
Sample usage:\n
python extractlikwidmetrics.py -path \'examples/151217_phi1_node/' -prefix "EulerFlow"
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-path',required=True,help="Directory containing the Peano output files.")
parser.add_argument('-prefix',required=True,help="Prefix of the Peano output files.")

args     = parser.parse_args();

root_dir = args.path
prefix   = args.prefix

extract_likwid_metrics(root_dir,prefix)
sort_table(prefix+".csv")
