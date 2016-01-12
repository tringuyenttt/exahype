import argparse
from argparse import RawTextHelpFormatter

import matplotlib.pylab
import matplotlib.pyplot as plt

import runtimeparser as rp
import hpclib as hpc
from plotting import scalingplot as sp

"""
.. module:: usertimeplot
  :platform: Unix, Windows
  :synopsis: Creates a speedup plot based on  Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Creates a speedup plot based on  Peano output files with specific file naming pattern.
"""

def plot_multithreading_adapter_scaling(root_dir,prefix,legend,adapters,process_counts,thread_counts,n_runs,cc,mode,ylim,per_iteration=False,hyperthreading=False,annotate=False):
    """
    Creates a scaling plot for the cumulative user time spent within the 
    specified adapters.
   
    Args:
      root_dir (str[]):
         Directories containing the Peano output files. (Implementation does currently only support one root dir element.)
      prefix (str[]):
         Prefix of the files - usually the date of the test and an identifier for the machine and the MPI process that has written the output files. Must be supplied once per 'root_dir' entry.
      legend (str[]):
         Legend entry for the each data set - usually an identifier for the machine and the MPI process that has written the output files. Must be supplied once per 'root_dir' entry.
      adapters (str[]):
         Name of the adapters. (Use 'Total' for the  cumulative time forall adapters. Does not make sense with per_iteration switched on.)
      process_counts (int[]):
         MPI process counts.
      thread_counts  (int[]):
         Threads per MPI process.
      n_runs (int):
         Number of runs for each 'n' and 't' combination [default=1].
      cc (str): 
         Compiler.
      mode (str):
         Shared memory mode.
      ylim (float):
         Upper limit for the y-Axis.
      hyperthreading (bool):
         The last thread count corresponds to a hyperthreading run.
      annotate (bool):
         Annotate the plots with the speedup values.
      per_iteration (bool):
         Use the adapter times per iteration.
    """
    n_process_counts = len(process_counts)
    n_thread_counts  = len(thread_counts)
    
    times                = rp.parse_all_adapter_times(root_dir[0],prefix[0],process_counts,thread_counts,n_runs,cc[0],mode[0],per_iteration)
    total_times        = rp.sum_all_adapter_times(times,n_process_counts,n_thread_counts)
    times['Total']    = total_times
    cumulative_times = rp.sum_adapter_times(times,adapters,n_process_counts,n_thread_counts)
    
    # Plotting
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    
    speedup_measured = hpc.compute_speedup_2(cumulative_times['avg_usertime'],cumulative_times['avg_usertime'][0][0])
    speedup_measured = speedup_measured[0];
    
    p         = map(int,thread_counts)
    p2        = p
    p_ticks = map(str,p)
    
    # Ideal speedup
    # In the hyper-threading case, we only want to plot a line
    # for the 'real' threads.
    if not hyperthreading: 
        plt.plot(p,p,label=r"ideal",markersize=4,marker="",markevery=1,lw=1.2,linestyle="dashed",color="grey")
    else:
        p2        = p[0:-1] # [inclusive:exclusive]
        p2.append(p2[-1]+2)
        p_ticks[-1] = "%d+HT" % p[-2]
        plt.plot(p[0:-1],p[0:-1],label=r"ideal",markersize=4,marker="",markevery=1,lw=1.2,linestyle="dashed",color="grey")
    
    # Measured speedup
    sp.plot_scaling(ax,p,speedup_measured,legend[0],"blue","s",hyperthreading,annotate)

    plt.ylabel(r"speedup", fontsize=12)    
    plt.xlabel(r"number of cores", fontsize=12)
    plt.grid(True)

    ax.set_xlim(0.8,p2[-1]+0.2)
    plt.xticks(p2,p_ticks)
    plt.yticks(p,map(str,p))
    ax.set_ylim(0.8,ylim+0.2)

    plt.suptitle("",fontsize=12)
    plt.tick_params(axis="both", which="major", labelsize=10)
    plt.tick_params(axis="both", which="minor", labelsize=10)
    
    fig = plt.gcf()
    matplotlib.pylab.legend(loc="best",fontsize="10")	
    DefaultSize = matplotlib.pylab.gcf().get_size_inches()
    fig.set_size_inches( (DefaultSize[0]/10, DefaultSize[1]/10) )
    fig.set_size_inches(7.25,7.25)
    plt.savefig("%s/%s_%s.pdf" % (root_dir[0],prefix[0],'+'.join(adapters)), bbox_inches="tight")
    plt.show()
    return

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'plot_multithreading_adapter_scaling' above.
help = '''
Creates a speedup plot based on Peano output files with specific file naming pattern.
If multiple adapters are specified, then the cumulative user times are used to compute the speedup.
\n\n
Sample usage:\n
python usertimeplot.py -path \'examples/151217_phi1_node\' -prefix \'151217\' -legend \'2x Xeon  5-2650 @ 2.00GHz\' -mode tbb -cc icpc -ylim 16 -per_iteration -adapter \'Predictor\' \'Corrector\' -t 1 2 4 6 8 10 12 16 32 -hyperthreading'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument("-path",nargs='+',required=True,help='Directories containing the Peano output files. (Implementation does currently only support one directory.)')
parser.add_argument("-prefix",default=[''],nargs='+',required=True,help='Prefix of the files - usually the date of the test and an identifier for the machine and the MPI process that has written the output files. Must be supplied once per \'root_dir\' entry.')
parser.add_argument("-legend",nargs='+',required=True,help='Legend entry for the each data set - usually an identifier for the machine and the MPI process that has written the output files. Must be supplied once per \'root_dir\' entry.')
parser.add_argument("-adapter",default=['Total'],nargs='+',help='Name of the adapters. (Use \'Total\' for the  cumulative time forall adapters. Does not make sense with per_iteration switched on.)')
parser.add_argument("-n",default=[1],nargs='+',help='MPI process counts [default=1].')
parser.add_argument("-t",default=[1],nargs='+',required=True,help='Threads per MPI process [default=1].')
parser.add_argument("-r",default=1,help='Number of runs for all \'n\' and \'t\' combinations [default=1].')
parser.add_argument("-cc",default='icpc',nargs='+',help='Compiler [default=\'icpc\']')
parser.add_argument("-mode",default='tbb',nargs='+',help='Shared memory mode [default=\'tbb\']')
parser.add_argument("-ylim",required=True,help='Upper limit for the y-Axis.')
parser.add_argument('-hyperthreading', action='store_true', default=False,help='The last thread count corresponds to a hyperthreading run.')
parser.add_argument('-annotate', action='store_true', default=False,help='Annotate the plots with the speedup values.')
parser.add_argument('-per_iteration', action='store_true', default=False,help='Use the adapter times per iteration instead of the total times.')

args           = parser.parse_args();
root_dir       = args.path
prefix         = args.prefix
legend         = args.legend
adapter        = args.adapter
process_counts = args.n
thread_counts  = args.t
n_runs         = args.r
cc             = args.cc
mode           = args.mode
ylim           = float(args.ylim)
hyperthreading = args.hyperthreading
annotate       = args.annotate
per_iteration  = args.per_iteration

plot_multithreading_adapter_scaling(root_dir,prefix,legend,adapter,process_counts,thread_counts,n_runs,cc,mode,ylim,per_iteration,hyperthreading,annotate)
