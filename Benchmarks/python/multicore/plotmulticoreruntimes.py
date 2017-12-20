#!/usr/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import matplotlib.pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter

# FROM: https://stackoverflow.com/questions/17687213/how-to-obtain-the-same-font-style-size-etc-in-matplotlib-output-as-in-latex
#Direct input 
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params) 

import operator
import csv


def column(matrix, i):
    return [row[i] for row in matrix]

def query_table(table,mesh,order,algorithm,cores): 
    '''
    Read a table and filter out certain
    columns 
    Args:
    table
       name of the csv file
    mesh
       mesh identifier; something like 'regular-0', 'regular-1'
    order
       approximaation order
    algorithm
       'fused' or 'nonfused'
    cores
       list of cores
    '''
    # print(table); print(mesh); print(order); print(algorithm); print(cores)
    dataFile = open(table, 'r')
    reader   = csv.reader(dataFile,delimiter='&')
    result   = list(filter(lambda x : x[0]==mesh and int(x[1])==order and x[4]==algorithm and int(x[7]) in cores, reader))
    dataFile.close() 
    
    # result row 10 (starting with 0) contains the total runtime (per iteration)
    return result


'''
.. module:: usertimeplot
  :platform: Unix, Windows, Mac
  :synopsis: Creates a speedup plot based on Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Creates a speedup plot based on user and CPU times stored in CSV tables.
'''

def plot_time_per_timestep(table,mesh,algorithms,cores,orders,y,ymin,ymax):
    '''
    '''
    # Plot
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    
    colors           = ['k','0.33']    # per algorithm
    markerFaceColors = ['k','None']       # per algorithm
    markers          = ['o','^','s','D']  # per order
    
    for ialg in range(0,len(algorithms)):
        algorithm = algorithms[ialg]
        markerSpacing=3
        if ialg==1:
            markerSpacing=4
        
        for ip in range(0,len(orders)):
            p = orders[ip]
            identifier=r'$p=%d\,\rm{(%s)}$' % (p,algorithm)
            if ialg*len(orders)+ip > len(orders):
                identifier=None
            plt.plot(cores,y[algorithm][ip],lw=1.1,markersize=5,color=colors[ialg],markeredgecolor=colors[ialg],markerfacecolor=markerFaceColors[ialg],marker=markers[ip],label=identifier,markevery=markerSpacing)
        
    plt.ylabel(r'time per cell per real. timestep [s]', fontsize=8)
    plt.xlabel(r'cores', fontsize=8)
    plt.grid(True, which='both')

    xticks = cores
    plt.tick_params(axis='both', which='major', labelsize=7)
    plt.axes().xaxis.set_major_locator(FixedLocator(cores[0::2]))
    plt.axes().xaxis.set_minor_locator(FixedLocator(cores[1::2]))
    
    ax.yaxis.grid(False, which='minor')
    ax.set_ylim(ymin*0.9,ymax*1.1)
    
    # Write files
    fig = plt.gcf()
    plt.legend(loc='upper right',fontsize=7,framealpha=0.5)
    fig.set_size_inches(2.40,2.20) # width: 0.470 * SIAM SISC \textwidth (=5.125in)
    
    # linear
    filename = table.replace('.runtimes.csv','-runtimes-'+mesh+'-linear')
    ax.set_yscale('linear')
    plt.savefig('%s.pdf' % filename, bbox_inches='tight')
    plt.savefig('%s.png' % filename, bbox_inches='tight')
    # log2
    filename = table.replace('.runtimes.csv','-runtimes-'+mesh+'-log2')
    ax.set_yscale('log',basey=2)
    plt.savefig('%s.pdf' % filename, bbox_inches='tight')
    plt.savefig('%s.png' % filename, bbox_inches='tight')
    # log10
    filename = table.replace('.runtimes.csv','-runtimes-'+mesh+'-log10')
    ax.set_yscale('log',basey=10)
    plt.savefig('%s.pdf' % filename, bbox_inches='tight')
    plt.savefig('%s.png' % filename, bbox_inches='tight')
    print('PDF and PNG output written.')

def plot_updates_per_timestep(table,mesh,algorithms,cores,orders,y,ymin,ymax):
    '''
    '''
    y_inv = y
    for ip in range(0,len(orders)):
        for algorithm in algorithms:
            for i in range(0,len(y[algorithm][ip])):
               y_inv[algorithm][ip][i]  =  1.0 / y[algorithm][ip][i] 
    
    # Plot
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    
    colors           = ['k','0.33']    # per algorithm
    markerFaceColors = ['k','None']       # per algorithm
    markers          = ['o','^','s','D']  # per order
    
    for ialg in range(0,len(algorithms)):
        algorithm = algorithms[ialg]
        markerSpacing=3
        if ialg==1:
            markerSpacing=4
        
        for ip in range(0,len(orders)):
            p = orders[ip]
            identifier=r'$p=%d\,\rm{(%s)}$' % (p,algorithm)
            if ialg*len(orders)+ip > len(orders):
                identifier=None
            plt.plot(cores,y_inv[algorithm][ip],lw=1.1,markersize=5,color=colors[ialg],markeredgecolor=colors[ialg],markerfacecolor=markerFaceColors[ialg],marker=markers[ip],label=identifier,markevery=markerSpacing)
    
    plt.ylabel(r'cell updates per real. timestep [s]', fontsize=8)
    plt.xlabel(r'cores', fontsize=8)
    plt.grid(True, which='both')
    
    xticks     = cores
    plt.tick_params(axis='both', which='major', labelsize=7)
    plt.axes().xaxis.set_major_locator(FixedLocator(cores[0::2]))
    plt.axes().xaxis.set_minor_locator(FixedLocator(cores[1::2]))
    
    ax.yaxis.grid(False, which='minor')
    ax.set_ylim(1.0/ymax*0.9,1.0/ymin*1.1)
    
    # Write files
    fig = plt.gcf()
    plt.legend(loc='lower right',fontsize=7,framealpha=0.5)
    fig.set_size_inches(2.40,2.20) # width: 0.470 * SIAM SISC \textwidth (=5.125in)
    
    # linear
    filename = table.replace('.runtimes.csv','-updates-'+mesh+'-linear')
    ax.set_yscale('linear')
    plt.savefig('%s.pdf' % filename, bbox_inches='tight')
    plt.savefig('%s.png' % filename, bbox_inches='tight')
    # log2
    filename = table.replace('.runtimes.csv','-updates-'+mesh+'-log2')
    ax.set_yscale('log',basey=2)
    plt.savefig('%s.pdf' % filename, bbox_inches='tight')
    plt.savefig('%s.png' % filename, bbox_inches='tight')
    # log10
    filename = table.replace('.runtimes.csv','-updates-'+mesh+'-log10')
    ax.set_yscale('log',basey=10)
    plt.savefig('%s.pdf' % filename, bbox_inches='tight')
    plt.savefig('%s.png' % filename, bbox_inches='tight')
    print('PDF and PNG output written.')

def plot_runtimes(table,mesh,orders,dim,ignoreReruns):
    '''
    Creates a scaling plot for the cumulative user time spent within the 
    specified adapters.
   
    Args:
      table (str):
         A csv file
    '''
    mesh2cells = {}
    mesh2cells['regular-0'] = 27**dim
    mesh2cells['regular-1'] = 81**dim
    mesh2cells['regular-2'] = 243**dim
    
    algorithms = ['fused', 'nonfused']
    #cores      = [1, 2, 3, 4, 5, 6 ,7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 ,24]
    cores      = [1, 2, 3, 4, 6 ,8, 9, 12, 15, 16, 18, 21, 24]
    nOrders    = len(orders)
    nCores     = len(cores)
    
    # Read data from table and process
    ymin = 10**20
    ymax = 0
    
    y = {}
    y['fused']   ={}
    y['nonfused']={}
    for ip in range(0,len(orders)):
        p = orders[ip]
        
        dataFused    = query_table(table,mesh,p,'fused',cores)
        dataNonfused = query_table(table,mesh,p,'nonfused',cores)
        
        y['nonfused'][ip]  = list( map( float, column(dataNonfused, 10) ) )
        if ignoreReruns is True:
            y['fused'][ip] = list( map( float, column(dataFused, 11) ) )
        else:
            y['fused'][ip] = list( map( float, column(dataFused, 10) ) )
        # Normalise fused values w.r.t. number of cells
        itFused    = float(dataFused[0][9])
        itNonfused = float(dataNonfused[0][9])
        for i in range(0,len(y['fused'][ip])):
            y['fused'][ip][i] /= mesh2cells[mesh]
            y['fused'][ip][i] *= itFused / itNonfused # Normalise fused values wrt iterations of nonfused
        for i in range(0,len(y['nonfused'][ip])):
            y['nonfused'][ip][i] /=mesh2cells[mesh]
        
        for algorithm in algorithms:
            ymax = max( ymax, max(y[algorithm][ip]) )
            ymin = min( ymin, min(y[algorithm][ip]) )
    
    # Plot    
    plot_time_per_timestep(table,mesh,algorithms,cores,orders,y,ymin,ymax)
    plot_updates_per_timestep(table,mesh,algorithms,cores,orders,y,ymin,ymax)

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'plot_multithreading_adapter_scaling' above.
help = '''
Creates a speedup plot based on the given tables containing adapter user and cpu times.
If multiple adapters are specified, then the cumulative user times are used to compute the speedup.
\n\n
Sample usage:\n
python plotmulticorespeedup.py -table 
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-table',required=True,help='A .runtimes.csv file created with the extractruntimes.py script.')
#parser.add_argument('-prefix',required=True,help='Prefix of the plot files.')
parser.add_argument('-mesh',required=True,help='Mesh identifier.')
parser.add_argument('-p',nargs='+',required=True,help='Polynomial orders (multiple values possible).')
parser.add_argument('-dim',required=True,help='Space dimensions.')
parser.add_argument('--ignore-reruns',dest='ignoreReruns',required=False,action='store_true',help="Do not consider predictor reruns in the runtime calculations.")
parser.set_defaults(ignoreReruns=False)

args           = parser.parse_args();

table          = args.table
mesh           = args.mesh
p              = args.p
p              = list(map(int,p))
dim            = int(args.dim)
ignoreReruns   = args.ignoreReruns


plot_runtimes(table,mesh,p,dim,ignoreReruns)
