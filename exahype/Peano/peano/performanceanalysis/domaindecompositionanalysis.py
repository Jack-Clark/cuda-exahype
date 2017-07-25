import argparse
from argparse import RawTextHelpFormatter

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!

import pylab
import networkx
import datetime
import performanceanalysisroutines


########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documen\tation
# of the function 'plot_multithreading_adapter_scaling' above.
help = '''
Creates a brief overview of the domain decomposition
\n\n
Sample usages:\n
python domain-decomposition-analysis.py -file 112x16-0.results -dimension 2 -domainoffset 0.0 0.0 -domainsize 15.0 15.0\n
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-file',required=True,help="Input file")
parser.add_argument('-dimension',required=True,help="Dimension of problem. Either 2 or 3.")
parser.add_argument('-domainoffset',nargs="+",required=True,help="Offset of bounding box.")
parser.add_argument('-domainsize',nargs="+",required=True,help="Size of domain's bounding box.")
args   = parser.parse_args();


outputFileName  = args.file + ".html"
outFile        = open( outputFileName, "w" )
outFile.write( 
  "<html>\
     <body>\
     <h1>Peano Domain Decomposition Report</h1>\
     <p>Report time stamp: " + str(datetime.datetime.utcnow()) + "</p>\
     <p>Author: Tobias Weinzierl</p>\
     <p>Data file: " + args.file + "</p>\
  ")


dim = int(args.dimension)


numberOfRanks = performanceanalysisroutines.getNumberOfRanks(args.file)


(parents,levels,offset,volume) = performanceanalysisroutines.plotLogicalTopology(args.file,numberOfRanks,dim);
outFile.write( "<h2>Logical topology</h2>" )
outFile.write( "<a href=\"" + args.file + ".topology.large.png\"><img src=\"" + args.file + ".topology.png\" /></a>" )


(volumes,overlaps,work) = performanceanalysisroutines.computeVolumesOverlapsWork(numberOfRanks,volume,offset,dim,args.domainoffset,args.domainsize,parents)



print "validate data ",
for i in range(0,numberOfRanks):
  print ".",
  if overlaps[i]<0.0:
    print "\nERROR for rank " + str(i) +  \
          ": region=" + str(offset[i]) + "x" + str(volume[i]) + \
          ", bounding box=" + str(args.domainoffset) + "x" + str(args.domainsize) + \
          ", volume= " + str(volumes[i]) + ", overlaps=" + str(overlaps[i]) 
    print "overlaps may not be negative"
  if volumes[i]<0.0:
    print "\nERROR for rank " + str(i) +  \
          ": region=" + str(offset[i]) + "x" + str(volume[i]) + \
          ", bounding box=" + str(args.domainoffset) + "x" + str(args.domainsize) + \
          ", volume= " + str(volumes[i]) + ", overlaps=" + str(overlaps[i]) 
    print "volumes may not be negative"
  if overlaps[i]>volume[i]:
    print "\nERROR for rank " + str(i) +  \
          ": region=" + str(offset[i]) + "x" + str(volume[i]) + \
          ", bounding box=" + str(args.domainoffset) + "x" + str(args.domainsize) + \
          ", volume= " + str(volumes[i]) + ", overlaps=" + str(overlaps[i]) 
    print "overlaps have to be bigger than volume"
  if work[i]<0:
    print "\nERROR for rank " + str(i) +  \
          ": region=" + str(offset[i]) + "x" + str(volume[i]) + \
          ", bounding box=" + str(args.domainoffset) + "x" + str(args.domainsize) + \
          ", volume= " + str(volumes[i]) + ", overlaps=" + str(overlaps[i]) + ", work=" + str(work[i])
    print "work has to be positive"
print " done "



performanceanalysisroutines.plotWorkloadAndResponsibilityDistribution(numberOfRanks,volumes,overlaps,work,args.file);



outFile.write( "<h2>Work statistics</h2>" )
outFile.write( "<a href=\"" + args.file + ".work-distribution.large.png\"><img src=\"" + args.file + ".work-distribution.png\" /></a>" )
outFile.write( "<a href=\"" + args.file + "-symlog.work-distribution.large.png\"><img src=\"" + args.file + "-symlog.work-distribution.png\" /></a>" )
outFile.write( "<p>The filled region is the actual local work volume of a rank. It has to be smaller than the region of responsibility that might overlap the actual domain.</p>" )


if dim==2:
 outFile.write( "<h2>Domain decomposition (level by level)</h2>" )
 for l in range(1,max(levels)+1):
  performanceanalysisroutines.plot2dDomainDecompositionOnLevel(l,numberOfRanks,args.domainoffset,args.domainsize,offset,volume,levels,args.file)
  outFile.write( "<a href=\"" + args.file + ".level" + str(l) + ".pdf\"> <img src=\"" + args.file + ".level" + str(l) + ".png\" /> </a> " )
 
 
outFile.write( "</html>" )

  