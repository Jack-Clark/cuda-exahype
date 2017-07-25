import sys
import re

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!

import pylab
import networkx 


def getNumberOfRanks(filename):
  print "parse input file header ",
  try:
    inputFile = open( filename,  "r" )
    for line in inputFile:
      print ".",
      if "processes:" in line:
        numberOfRanks = int(line.split( "processes:" )[1].split( "," )[0])
        print str(numberOfRanks) + " MPI ranks used"
        return numberOfRanks

  except Exception as inst:
    print "failed to read " + filename
    print inst
  return 0
  
  
def getBoundingBox(line):
  line    = line.replace( "]x[", "," )
  results = line.split( "node for subdomain [")[1].split( "] on level" )[0].split( "," )
  return results

  
searchPatternNumberOfCells         = "peano::performanceanalysis::DefaultAnalyser::endIteration.*cells=\("
searchPatternNumberOfLocalVertices = "peano::performanceanalysis::DefaultAnalyser::endIteration.*local-vertices="
searchPatternTTotal                = "([0-9]\.?[0-9]*).*peano::performanceanalysis::DefaultAnalyser::endIteration.*t_total=\("
searchPatternTTraversal            = "peano::performanceanalysis::DefaultAnalyser::endIteration.*t_traversal=\("



numberOfInnerLeafCells  = {}
numberOfInnerCells      = {}
numberOfOuterLeafCells  = {}
numberOfOuterCells      = {}
numberOfLocalCells      = {}
numberOfLocalVertices   = {}
tTotal                  = {}
tTraversal              = {}

  
  
#
# Display constants. Might be reset later
#
AlphaValue = 1.0


#
# Use these pairs also for sender and receiver
#
class Pair:
  def __init__(self, master, worker):
    self.master         = master
    self.worker         = worker
    self.count          = 0
    self.maxTime        = 0.0
    self.average        = 0
    self.maxCardinality = 0


def switchToLargePlot():
  DefaultSize = pylab.gcf().get_size_inches()
  pylab.gcf().set_size_inches( (DefaultSize[0]*10, DefaultSize[1]*10) )





def getNumberOfThreads(filename):
  print "parse input file header ",
  try:
    inputFile = open( filename,  "r" )
    for line in inputFile:
      print ".",
      if "threads:" in line:
        numberOfThreads = int(line.split( "threads:" )[1].split( "(per process)" )[0])
        print str(numberOfThreads) + " threads ranks used"
        return numberOfThreads

  except Exception as inst:
    print "failed to read " + filename
    print inst
  return 0



def switchBackToStandardPlot():
  DefaultSize = pylab.gcf().get_size_inches()
  pylab.gcf().set_size_inches( (DefaultSize[0]/10, DefaultSize[1]/10) )


def drawTreeGraph(myGraph):
  try:
    from networkx import graphviz_layout
    pos=networkx.graphviz_layout(myGraph,prog='twopi',args='')
  except:
    try:
      pos=networkx.graphviz_layout(myGraph,prog='dot',args='')
    except:
      pos=networkx.spring_layout(myGraph)
      print "fall back to spring layout. Nicer graphs might result if either compatible PyGraphviz or Pydot packages were available"
  networkx.draw(
    myGraph,
    pos,
    with_labels=True,
    arrows=True,
    node_color='#667766',
    node_size=10,
    alpha=0.2
  )


def parseInputFile(numberOfRanks,inputFileName):
  print "prepare data structures",
  
  myNumberOfRanks = numberOfRanks
  if myNumberOfRanks==0:
    myNumberOfRanks = 1
  for rank in range(0,myNumberOfRanks):
    numberOfInnerLeafCells[rank]  = []
    numberOfInnerCells[rank]      = []
    numberOfOuterLeafCells[rank]  = []
    numberOfOuterCells[rank]      = []
    numberOfLocalCells[rank]      = []
    numberOfLocalVertices[rank]   = []
    tTotal[rank]                  = []
    tTraversal[rank]              = []
    print ".",
  print " done"
   
  try:
    inputFile = open( inputFileName,  "r" )
    print "parse ",
    for line in inputFile:
        if ("DefaultAnalyser" in line):
          m = re.search( searchPatternNumberOfCells, line )
          if (m):
            rank  = 0
            if numberOfRanks>0:
              rank  = int(line.split( "rank:" )[-1].split( " " )[0])
            token = line.replace("(",",").replace(")","").strip().split(",")
            numberOfInnerLeafCells[rank].append(float(token[-5]))
            numberOfOuterLeafCells[rank].append(float(token[-4]))
            numberOfInnerCells[rank].append(float(token[-3]))
            numberOfOuterCells[rank].append(float(token[-2]))
            numberOfLocalCells[rank].append(float(token[-1]))
            print ".",
          m = re.search( searchPatternNumberOfLocalVertices, line )
          if (m):
            rank  = 0
            if numberOfRanks>0:
              rank  = int(line.split( "rank:" )[-1].split( " " )[0])
            token = line.strip().split("=")
            numberOfLocalVertices[rank].append(float(token[-1]))
            print ".",
          m = re.search( searchPatternTTotal, line )
          if (m):
            rank  = 0
            if numberOfRanks>0:
              rank  = int(line.split( "rank:" )[-1].split( " " )[0])
            timeStamp = float( m.group(1) )
            tTotal[rank].append(timeStamp)
            print ".",
          m = re.search( searchPatternTTraversal, line )
          if (m):
            rank  = 0
            if numberOfRanks>0:
              rank  = int(line.split( "rank:" )[-1].split( " " )[0])
            token = line.strip().split("=(")[-1].split(",")[0]
            tTraversal[rank].append(float(token))
            print ".",
          
    print " done"
  except Exception as inst:
    print "failed to read " + inputFileName
    print inst



def plotConcurrency(rank,inputFileName):
  pylab.clf()
  #pylab.title( "Concurrency levels" )

  ax = pylab.gca()
  #ax.set_aspect('equal','box')
  
  timeStampPattern = "([0-9]+\.?[0-9]*)"
  floatPattern = "([0-9]\.?[0-9]*)"
  
  searchPattern              = timeStampPattern + ".*rank:" + str(rank) + " .*peano::performanceanalysis::DefaultAnalyser::changeConcurrencyLevel" + \
                               ".*dt=" + floatPattern + \
                               ".*cpu-time=" + floatPattern + \
                               ".*concurrent-time=" + floatPattern + \
                               ".*potential-concurrent-time=" + floatPattern + \
                               ".*max-concurrency-level=" + floatPattern + \
                               ".*max-potential-concurrency-level=" + floatPattern + \
                               ".*background-tasks=" + floatPattern

  lastTimeStamp = 0
  
  
  maxPotentialConcurrencyLevelLastPoint = 0
  potentialConcurrencyTimeLastPoint     = 0
  realConcurrencyLastPoint              = 1
  
  try:
    inputFile = open( inputFileName,  "r" )
    print "parse concurrency level of rank " + str(rank)
    for line in inputFile:
      m = re.search( searchPattern, line )
      if (m):
        timeStamp                      = float( m.group(1) )
        dt                             = float( m.group(2) )
        cpuTime                        = float( m.group(3) ) 
        concurrentTime                 = float( m.group(4) )
        potentialConcurrentTime        = float( m.group(5) ) 
        maxConcurrencyLevel            = float( m.group(6) )
        maxPotentialConcurrencyLevel   = float( m.group(7) ) 
        backgroundTasks                = float( m.group(8) )
        
        
        #print "found dt=" + str(dt) + ", cpuTime=" + str(cpuTime)
         
  
        backgroundTasksBar         = pylab.Rectangle([lastTimeStamp,maxConcurrencyLevel],timeStamp-lastTimeStamp,backgroundTasks,facecolor="#0000ff",edgecolor="#0000ff")
        maxConcurrencyLevelBar     = pylab.Rectangle([lastTimeStamp,0],timeStamp-lastTimeStamp,maxConcurrencyLevel,              facecolor="#ff0000",edgecolor="#ff0000")
        concurrentTimeBar          = pylab.Rectangle([lastTimeStamp,0],timeStamp-lastTimeStamp,concurrentTime/dt,                facecolor="#aa0000",edgecolor="#bb3344")

        maxPotentialConcurrencyLevelSymbol = pylab.plot([lastTimeStamp,timeStamp],[maxPotentialConcurrencyLevelLastPoint, maxPotentialConcurrencyLevel], "-",  color="#00ff00" )
        potentialConcurrencyTimeSymbol     = pylab.plot([lastTimeStamp,timeStamp],[potentialConcurrencyTimeLastPoint,     potentialConcurrentTime/dt]  , ":",  color="#008800" )
        realConcurrencyLevelSymbol         = pylab.plot([lastTimeStamp,timeStamp],[realConcurrencyLastPoint,              cpuTime/dt                ]  , "-",  color="#000000")

        maxPotentialConcurrencyLevelLastPoint = maxPotentialConcurrencyLevel
        potentialConcurrencyTimeLastPoint     = potentialConcurrentTime/dt
        realConcurrencyLastPoint              = cpuTime/dt

        ax.add_patch(backgroundTasksBar)
        ax.add_patch(maxConcurrencyLevelBar)
        ax.add_patch(concurrentTimeBar)

        lastTimeStamp = timeStamp
        print ".",
    #pylab.plot([0,timeStamp],[1,1], "--", color="#000000")
  except Exception as inst:
    print "failed to read " + inputFileName
    print inst
  
  ax.autoscale_view()
  ax.set_yscale('symlog', basey=2)
  
  pylab.yticks( 
    [1,2,8,12,16,18,24,60,72,120,180,240,480], 
    ["1","2","8","12","16","18","24","60","72","120","180","240","480"] 
  )
  
  pylab.xlabel('t')
  pylab.ylabel('Concurrency level')
  pylab.grid(True)
  pylab.savefig( inputFileName + ".html-rank-" + str(rank) + ".concurrency.png" )
  pylab.savefig( inputFileName + ".html-rank-" + str(rank) + ".concurrency.pdf" )
  switchToLargePlot()
  pylab.savefig( inputFileName + ".html-rank-" + str(rank) + ".concurrency.large.png" )
  pylab.savefig( inputFileName + ".html-rank-" + str(rank) + ".concurrency.large.pdf" )
  switchBackToStandardPlot()  





#
# Creates the one big trace picture where we see what different ranks do at 
# different times
#
# @writes A file inputFileName.mpi-phases.large.png
# @writes A file inputFileName.mpi-phases.large.pdf
# @writes A file inputFileName.mpi-phases.png
# @writes A file inputFileName.mpi-phases.pdf
#
#
def plotMPIPhases(numberOfRanks,inputFileName):
  beforeInTraversalColor  = "#ff3434"
  inTraversalColor        = "#00ab00"
  afterInTraversalColor   = "#560000"
  afterBoundaryExchange   = "#0000ab"
  prepareAsynchronousBoundaryExchangeColor = "#ffff00"
  releaseAsynchronousBoundaryExchangeColor = "#abab00"

  pylab.clf()
  #pylab.figure(figsize=(numberOfRanks,numberOfRanks))
  pylab.title( "MPI phases overview" )
  ax = pylab.gca()
  #ax.set_aspect('equal','box')
  
  timeStampPattern = "([0-9]+\.?[0-9]*)"
  floatPattern = "([0-9]\.?[0-9]*)"
  
  beginIterationPattern      = timeStampPattern + ".*rank:(\d+)*.*peano::performanceanalysis::DefaultAnalyser::beginIteration"
  enterCentralElementPattern = timeStampPattern + ".*rank:(\d+)*.*peano::performanceanalysis::DefaultAnalyser::enterCentralElementOfEnclosingSpacetree"
  leaveCentralElementPattern = timeStampPattern + ".*rank:(\d+)*.*peano::performanceanalysis::DefaultAnalyser::leaveCentralElementOfEnclosingSpacetree.*t_central-tree-traversal=\(" + floatPattern
  endIterationPattern        = timeStampPattern + ".*rank:(\d+)*.*peano::performanceanalysis::DefaultAnalyser::endIteration.*t_traversal=\(" + floatPattern
  endDataExchange            = timeStampPattern + ".*rank:(\d+)*.*peano::performanceanalysis::DefaultAnalyser::endReleaseOfBoundaryData"
  prepareAsynchronousHeap    = timeStampPattern + ".*rank:(\d+)*.*peano::performanceanalysis::DefaultAnalyser::endToPrepareAsynchronousHeapDataExchange"
  releaseAsynchronousHeap    = timeStampPattern + ".*rank:(\d+)*.*peano::performanceanalysis::DefaultAnalyser::endToReleaseSynchronousHeapData"

  lastTimeStamp  = [0] * numberOfRanks

  def plotMPIPhasesBar(rank,timeStamp,color):
    if (lastTimeStamp[rank]==0):
      lastTimeStamp[rank] = timeStamp
    rectLength = timeStamp-lastTimeStamp[rank]
    rect = pylab.Rectangle([lastTimeStamp[rank],rank-0.5],rectLength,1,facecolor=color,edgecolor=color,alpha=Alpha)
    ax.add_patch(rect)
    lastTimeStamp[rank] = lastTimeStamp[rank] + rectLength
  
  Alpha = 0.5
  
  try:
    inputFile = open( inputFileName,  "r" )
    print "parse mpi phases",
    for line in inputFile:
      m = re.search( beginIterationPattern, line )
      if (m):
        rank = int( m.group(2) )
        timeStamp = float( m.group(1) )
        lastTimeStamp[rank] = timeStamp
        print ".",
        if (rank==0):
          pylab.plot((timeStamp, timeStamp), (-0.5, numberOfRanks+1), '--', color="#445544", alpha=Alpha)
        pylab.plot((timeStamp, timeStamp), (rank-0.5, rank+0.5), '-', color="#000000" )
      m = re.search( prepareAsynchronousHeap, line )
      if (m):
        rank = int( m.group(2) )
        timeStamp = float( m.group(1) )
        plotMPIPhasesBar(rank,timeStamp,prepareAsynchronousBoundaryExchangeColor)
      m = re.search( releaseAsynchronousHeap, line )
      if (m):
        rank = int( m.group(2) )
        timeStamp = float( m.group(1) )
        plotMPIPhasesBar(rank,timeStamp,releaseAsynchronousBoundaryExchangeColor)
      m = re.search( enterCentralElementPattern, line )
      if (m):
        rank = int( m.group(2) )
        timeStamp = float( m.group(1) )
        plotMPIPhasesBar(rank,timeStamp,beforeInTraversalColor)
      m = re.search( leaveCentralElementPattern, line )
      if (m):
        rank = int( m.group(2) )
        timeStamp = float( m.group(1) )
        plotMPIPhasesBar(rank,timeStamp,inTraversalColor)
      m = re.search( endIterationPattern, line )
      if (m):
        rank = int( m.group(2) )
        timeStamp = float( m.group(1) )
        plotMPIPhasesBar(rank,timeStamp,afterInTraversalColor)
      m = re.search( endDataExchange, line )
      if (m):
        rank = int( m.group(2) )
        timeStamp = float( m.group(1) )
        plotMPIPhasesBar(rank,timeStamp,afterBoundaryExchange)

    print " done"
  except Exception as inst:
    print "failed to read " + inputFileName
    print inst
  
  ax.invert_yaxis()
  ax.autoscale_view()
  pylab.xlabel('t')
  pylab.grid(False)
  pylab.savefig( inputFileName + ".mpi-phases.png" )
  pylab.savefig( inputFileName + ".mpi-phases.pdf" )
  switchToLargePlot()
  if numberOfRanks<=16:
    pylab.yticks([i for i in range(0,numberOfRanks)]) 
  else:
    pylab.yticks([i*16 for i in range(0,numberOfRanks/16)]) 
  pylab.savefig( inputFileName + ".mpi-phases.large.png" )
  pylab.savefig( inputFileName + ".mpi-phases.large.pdf" )
  switchBackToStandardPlot()  



def plotForkJoinStatistics(numberOfRanks,inputFileName):
  numberOfWorkingNodes = []
  numberOfIdleNodes    = []
  numberOfForks        = []
  numberOfJoins        = []
  timelineOfWorkingNodes = []
  timelineOfIdleNodes    = []
  timelineOfForks        = []
  timelineOfJoins        = []
  
  numberOfWorkingNodes.append(1)
  numberOfIdleNodes.append(numberOfRanks-1)
  numberOfForks.append(0)
  numberOfJoins.append(0)
  timelineOfWorkingNodes.append(0.0)
  timelineOfIdleNodes.append(0.0)
  timelineOfForks.append(0.0)
  timelineOfJoins.append(0.0)
  
  timeStampPattern  = "([0-9]+\.?[0-9]*)"
  
  searchPatternFork    = timeStampPattern + ".*peano::performanceanalysis::DefaultAnalyser::addWorker"
  searchPatternJoin    = timeStampPattern + ".*peano::performanceanalysis::DefaultAnalyser::removeWorker"
  
  try:
    inputFile = open( inputFileName,  "r" )
    print "parse fork/join statistics",
    for line in inputFile:
      m = re.search( searchPatternFork, line )
      if (m):
        timeStamp = float( m.group(1) )
        numberOfWorkingNodes.append(numberOfWorkingNodes[-1]+1)
        numberOfIdleNodes.append(numberOfIdleNodes[-1]-1)
        numberOfForks.append(numberOfForks[-1]+1)

        timelineOfWorkingNodes.append(timeStamp)
        timelineOfIdleNodes.append(timeStamp)
        timelineOfForks.append(timeStamp)
      m = re.search( searchPatternJoin, line )
      if (m):
        timeStamp = float( m.group(1) )
        
        numberOfWorkingNodes.append(numberOfWorkingNodes[-1]-1)
        numberOfIdleNodes.append(numberOfIdleNodes[-1]+1)
        numberOfJoins.append(numberOfJoins[-1]+1)

        timelineOfWorkingNodes.append(timeStamp)
        timelineOfIdleNodes.append(timeStamp)
        timelineOfJoins.append(timeStamp)
    print " done"
  except Exception as inst:
    print "failed to read " + inputFileName
    print inst
  
  
  pylab.clf()
  pylab.title( "Fork and join statistics" )

  pylab.plot(timelineOfWorkingNodes, numberOfWorkingNodes, 'o-',  markersize=10, color='#ffaa00', label='working nodes' )
  pylab.plot(timelineOfIdleNodes, numberOfIdleNodes,    '+-',  markersize=10, color='#00ffaa', label='idle nodes' )
  pylab.plot(timelineOfForks, numberOfForks,        '.-',  markersize=10, color='#aa00ff', label='total forks' )
  pylab.plot(timelineOfJoins, numberOfJoins,        'x-',  markersize=10, color='#ff00aa', label='total joins' )
  
  setGeneralPlotSettings()
  
  pylab.savefig( inputFileName + ".fork-join-statistics.png" )
  pylab.savefig( inputFileName + ".fork-join-statistics.pdf" )




  

def plotBoundaryLateSends(numberOfRanks,inputFileName):
  pairs = dict()
  totalMaxCardinality     = 0
  totalCount              = 0
  totalAverage            = 0
  waitingForNeighborLine  = "rank:(\d+)*.*peano::performanceanalysis::DefaultAnalyser::dataWasNotReceivedInBackground.*rank had to wait for (\d+) record\(s\)" \
                          + " from (\d+) on tag (\d) with page size (\d+)"

  print "search for late boundary answers ...",
  try:
    inputFile = open( inputFileName,  "r" )
    for line in inputFile:
      m = re.search( waitingForNeighborLine, line )
      if (m):
        master = int(m.group(1))
        worker = int(m.group(3))
        key = (master, worker)
        
        if not pairs.has_key(key):
          pairs[key] = Pair(master, worker)
        pair = pairs[key]
        
        cardinality = int(m.group(1))
        pair.count += 1
        totalCount += 1 
     
        pair.maxCardinality = max(pair.maxCardinality, cardinality)
        totalMaxCardinality = max(totalMaxCardinality, cardinality)

        pair.average += cardinality
        totalAverage += cardinality

    print "boundary analysis found " + str( len(pairs) ) + " entries of interest in trace file" 
  except Exception as inst:
    print "failed to read " + inputFileName
    print inst
  
  
  graph              = networkx.DiGraph()
  sparseAverageGraph = networkx.DiGraph()
  sparseMaxGraph     = networkx.DiGraph()
  for rank in range(0,numberOfRanks):
    graph.add_node(str(rank))
    sparseAverageGraph.add_node(str(rank))
    sparseMaxGraph.add_node(str(rank))
  for master in range(0,numberOfRanks):
    for worker in range(0,numberOfRanks):
      key = (master, worker)
      if pairs.has_key(key):
        pair = pairs[key]
        if pair.count>2:
          graph.add_edge(str(worker),str(master))
          #edge = pydot.Edge(str(sender),str(receiver), label="(" + str(pair.count) + "," + str(pair.maxCardinality) + "," + str(float(pair.average) / float(pair.count)) + ")", fontsize=str(myfontsize), labelfontcolor="blue" )
        
        if totalCount>0 and pair.count>2 and (float(pair.average) / float(pair.count) > float(totalAverage) / float(totalCount)):
          #edge = pydot.Edge(str(sender),str(receiver), label="(" + str(pair.count) + "," + str(pair.maxCardinality) + "," + str(float(pair.average) / float(pair.count)) + ")", fontsize=str(myfontsize), labelfontcolor="blue" )
          sparseAverageGraph.add_edge(str(worker),str(master))
 
        if pair.count>2 and ( float(pair.maxCardinality) > 0.9 * float(totalMaxCardinality) ):
          #edge = pydot.Edge(str(sender),str(receiver), label="(" + str(pair.count) + "," + str(pair.maxCardinality) + "," + str(float(pair.average) / float(pair.count)) + ")", fontsize=str(myfontsize), labelfontcolor="blue" )
          sparseMaxGraph.add_edge(str(worker),str(master))
  pylab.clf()
  pylab.title( "Late sends from neighbours" )
  drawTreeGraph(graph)
  pylab.savefig( inputFileName + ".boundary-data-exchange.png" )
  pylab.savefig( inputFileName + ".boundary-data-exchange.pdf" )
  switchToLargePlot()
  pylab.savefig( inputFileName + ".boundary-data-exchange.large.png" )
  pylab.savefig( inputFileName + ".boundary-data-exchange.large.pdf" )
  switchBackToStandardPlot()

  pylab.clf()
  pylab.title( "Late sends from neighbours (only edges with weight more than average)" )
  drawTreeGraph(sparseAverageGraph)
  pylab.savefig( inputFileName + ".boundary-data-exchange.sparse-average.png" )
  pylab.savefig( inputFileName + ".boundary-data-exchange.sparse-average.pdf" )
  switchToLargePlot()
  pylab.savefig( inputFileName + ".boundary-data-exchange.sparse-average.large.png" )
  pylab.savefig( inputFileName + ".boundary-data-exchange.sparse-average.large.pdf" )
  switchBackToStandardPlot()

  pylab.clf()
  pylab.title( "Late sends from neighbours (only 10% heaviest edges)" )
  drawTreeGraph(sparseMaxGraph)
  pylab.savefig( inputFileName + ".boundary-data-exchange.sparse-max.png" )
  pylab.savefig( inputFileName + ".boundary-data-exchange.sparse-max.pdf" )
  switchToLargePlot()
  pylab.savefig( inputFileName + ".boundary-data-exchange.sparse-max.large.png" )
  pylab.savefig( inputFileName + ".boundary-data-exchange.sparse-max.large.pdf" )
  switchBackToStandardPlot()


def plotMasterWorkerLateSends(numberOfRanks,inputFileName):
  GlobalSynchronisationOnRank0 = False
  
  #  <a href=\"" + outputFileName + ".master-worker-data-exchange.large.png\" /><img src=\"" + outputFileName + ".master-worker-data-exchange.png\" />\
  pairs = dict()
  totalMaxTime        = 0.0
  totalAverage        = 0
  totalCount          = 0
  floatPattern = "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
  waitingForWorkerLine = "rank:(\d+)*.*peano::performanceanalysis::DefaultAnalyser::endToReceiveDataFromWorker.+rank had to wait for worker (\d+) for ("+ floatPattern + ")s"

  print "search for late worker answers ",
  try:
    inputFile = open( inputFileName,  "r" )
    print "parse fork/join statistics",
    for line in inputFile:
      m = re.search( waitingForWorkerLine, line )
      if (m):
        master = int(m.group(1))
        worker = int(m.group(2))
        key = (master, worker)
        
        if not pairs.has_key(key):
          pairs[key] = Pair(master, worker)
        pair = pairs[key]
        
        time = float(m.group(3))
     
        pair.maxTime = max(pair.maxTime, time)
        totalMaxTime = max(totalMaxTime, time)
        pair.average += time
        totalAverage += time
        
        pair.count += 1
        totalCount += 1 
    print " found " + str( len(pairs) ) + " entries of interest in trace file" 
  except Exception as inst:
    print "failed to read " + inputFileName
    print inst
  
  
  graph              = networkx.DiGraph()
  sparseAverageGraph = networkx.DiGraph()
  sparseMaxGraph     = networkx.DiGraph()
  for rank in range(0,numberOfRanks):
    graph.add_node(str(rank))
    sparseAverageGraph.add_node(str(rank))
    sparseMaxGraph.add_node(str(rank))
  for master in range(0,numberOfRanks):
    for worker in range(0,numberOfRanks):
      key = (master, worker)
      if pairs.has_key(key):
        pair = pairs[key]
        if pair.count>2:
          graph.add_edge(str(worker),str(master))
          #edge = pydot.Edge(str(worker),str(master), label="(" + str(pair.count) + "," + str(pair.maxTime) + "," + str(float(pair.average) / float(pair.count)) + ")", fontsize=str(myfontsize), labelfontcolor="blue" )
        
        if totalCount > 0 and pair.count>2 and (float(pair.average) / float(pair.count) > float(totalAverage) / float(totalCount)):
          sparseAverageGraph.add_edge(str(worker),str(master))
          if master==0:
            GlobalSynchronisationOnRank0 = True
            print "detected very strong MPI synchronisation"
          #edge = pydot.Edge(str(worker),str(master), label="(" + str(pair.count) + "," + str(pair.maxTime) + "," + str(float(pair.average) / float(pair.count)) + ")", fontsize=str(myfontsize), labelfontcolor="blue" )
 
        if pair.count>2 and ( float(pair.maxTime) > 0.9 * float(totalMaxTime) ):
          sparseMaxGraph.add_edge(str(worker),str(master))
          if master==0:
            print "detected very strong MPI synchronisation"
            GlobalSynchronisationOnRank0 = True
          #edge = pydot.Edge(str(worker),str(master), label="(" + str(pair.count) + "," + str(pair.maxTime) + "," + str(float(pair.average) / float(pair.count)) + ")", fontsize=str(myfontsize), labelfontcolor="blue" )
          #sparseMaxGraph.add_edge(edge)      
  pylab.clf()
  pylab.title( "Late workers" )
  drawTreeGraph(graph)
  pylab.savefig( inputFileName + ".master-worker-data-exchange.png" )
  pylab.savefig( inputFileName + ".master-worker-data-exchange.pdf" )
  switchToLargePlot()
  pylab.savefig( inputFileName + ".master-worker-data-exchange.large.png" )
  pylab.savefig( inputFileName + ".master-worker-data-exchange.large.pdf" )
  switchBackToStandardPlot()

  pylab.clf()
  pylab.title( "Late workers (only edges with weight more than average)" )
  drawTreeGraph(sparseAverageGraph)
  pylab.savefig( inputFileName + ".master-worker-data-exchange.sparse-average.png" )
  pylab.savefig( inputFileName + ".master-worker-data-exchange.sparse-average.pdf" )
  switchToLargePlot()
  pylab.savefig( inputFileName + ".master-worker-data-exchange.sparse-average.large.png" )
  pylab.savefig( inputFileName + ".master-worker-data-exchange.sparse-average.large.pdf" )
  switchBackToStandardPlot()

  pylab.clf()
  pylab.title( "Late workers (only 10% heaviest edges)" )
  drawTreeGraph(sparseMaxGraph)
  pylab.savefig( inputFileName + ".master-worker-data-exchange.sparse-max.png" )
  pylab.savefig( inputFileName + ".master-worker-data-exchange.sparse-max.pdf" )
  switchToLargePlot()
  pylab.savefig( inputFileName + ".master-worker-data-exchange.sparse-max.large.png" )
  pylab.savefig( inputFileName + ".master-worker-data-exchange.sparse-max.large.pdf" )
  switchBackToStandardPlot()

  return GlobalSynchronisationOnRank0

#
# Plot topology to file and return an array where we bookmark for each rank what
# its parent is. The second entry of the return tuple is an array of the 
# corresponding levels.
#
def plotLogicalTopology(inputFileName,numberOfRanks,dim):
  parents  = [-1 for x in range(0,numberOfRanks)]
  levels   = [-1 for x in range(0,numberOfRanks)]
  if dim==2:
    volume = [(0.0,0.0) for x in range(0,numberOfRanks)]
    offset = [(0.0,0.0) for x in range(0,numberOfRanks)]
  else:
    volume = [(0.0,0.0,0.0) for x in range(0,numberOfRanks)]
    offset = [(0.0,0.0,0.0) for x in range(0,numberOfRanks)]
  level = [ 0 for x in range(0,numberOfRanks)]
  
  topologyGraph = networkx.MultiDiGraph()
  for rank in range(0,numberOfRanks):
    topologyGraph.add_node(str(rank))
    #topologyGraph.node[rank]['data']  = True
    #topologyGraph.node[rank]['label']  = str(rank)
  try:
    inputFile = open( inputFileName,  "r" )
    print "parse topology",
    for line in inputFile:
      if "start node for subdomain" in line:
        print ".",
        child = int(line.split( "rank:" )[1].split( " " )[0])
        parent = int(line.split( "with master" )[1])
        topologyGraph.add_edge(child,parent)
        parents[child]=parent
        fragments = getBoundingBox(line) 
        if dim==2:
          offset[child] = ( float(fragments[0]),float(fragments[1]) )
          volume[child] = ( float(fragments[2]),float(fragments[3]) )
        else:
          offset[child] = ( float(fragments[0]),float(fragments[1]),float(fragments[2]) )
          volume[child] = ( float(fragments[3]),float(fragments[4]),float(fragments[5]) )
        levels[child] = int(line.split( "on level" )[1].split( " with" )[0])
    print " done"
  except Exception as inst:
    print "failed to read " + inputFileName
    print inst
  pylab.clf()
  pylab.title( "Logical topology" )
  drawTreeGraph(topologyGraph)
  pylab.savefig( inputFileName + ".topology.png" )
  pylab.savefig( inputFileName + ".topology.pdf" )
  switchToLargePlot()
  pylab.savefig( inputFileName + ".topology.large.png" )
  pylab.savefig( inputFileName + ".topology.large.pdf" )
  switchBackToStandardPlot()
  return (parents,levels,offset,volume)


def setGeneralPlotSettings():
  pylab.xlabel('t')
  pylab.grid(True)
  # pylab.style.use('greyscale')
  # bmh              Hintergrund ist grau. Passt nicht so richtig hier. Im Gegensatz zu ggplot sind Gitterlinien gestrichelt
  # fivethirtyeight  Hintergrundfarbe ist homogen grau, kein Gitter. Es gibt einen grossen, dicken Rand. Passt zu Business-Grafiken 
  # ggplot           Hintergrund ist grau. Passt nicht so richtig hier 
  # greyscale        Geht bei mir nicht
  try:
    pylab.legend(fontsize=9, framealpha=0.5)
  except:
    # old pylab version
    l = pylab.legend(prop={'size':9})
    l.get_frame().set_alpha(0.5)



def plotGlobalGridOverview(numberOfRanks,inputFileName):
  pylab.clf()
  pylab.title( "Cells on global master" )

  pylab.plot(tTotal[0], numberOfInnerLeafCells[0], 'o-',  markersize=10, color='#ff0000', label='#inner leaf cells' )
  pylab.plot(tTotal[0], numberOfOuterLeafCells[0], '+-',  markersize=10, color='#00ff00', label='#outer leaf cells' )
  pylab.plot(tTotal[0], numberOfInnerCells[0],     '.-',  markersize=10, color='#aa0066', label='#inner cells' )
  pylab.plot(tTotal[0], numberOfOuterCells[0],     'x-',  markersize=10, color='#00aa66', label='#outer cells' )
  setGeneralPlotSettings()
  pylab.legend(fontsize=9, loc='upper left', framealpha=0.5)
  pylab.savefig( inputFileName + ".grid-overview-global-master.png" )
  pylab.savefig( inputFileName + ".grid-overview-global-master.pdf" )

  pylab.clf()
  pylab.title( "Local cells" )
  startRank = 1
  if (numberOfRanks==1):
    startRank = 0
  for rank in range(startRank,numberOfRanks):
    x = pylab.arange(0, len(numberOfLocalCells[rank]), 1.0)
    if len(tTotal[rank])==len(numberOfLocalCells[rank]):
      pylab.plot(tTotal[rank], numberOfLocalCells[rank], 'o',  color='#000000', alpha=AlphaValue, markersize=10)
  pylab.xlabel('t')
  pylab.savefig( inputFileName + ".local-cells.png" )
  pylab.savefig( inputFileName + ".local-cells.pdf" )

  globalNumberOfInnerLeafCells = [a for a in numberOfInnerLeafCells[0]]
  globalNumberOfOuterLeafCells = [a for a in numberOfOuterLeafCells[0]]
  globalNumberOfInnerCells     = [a for a in numberOfInnerCells[0]]
  globalNumberOfOuterCells     = [a for a in numberOfOuterCells[0]]
  
  for rank in range(1,numberOfRanks):
    currentEntry = len(globalNumberOfInnerLeafCells)-1
    for i in reversed(numberOfInnerLeafCells[rank]):
      globalNumberOfInnerLeafCells[currentEntry] = globalNumberOfInnerLeafCells[currentEntry] + i
      currentEntry = currentEntry - 1
    currentEntry = len(globalNumberOfOuterLeafCells)-1
    for i in reversed(numberOfOuterLeafCells[rank]):
      globalNumberOfOuterLeafCells[currentEntry] = globalNumberOfOuterLeafCells[currentEntry] + i
      currentEntry = currentEntry - 1
    currentEntry = len(globalNumberOfInnerCells)-1
    for i in reversed(numberOfInnerCells[rank]):
      globalNumberOfInnerCells[currentEntry] = globalNumberOfInnerCells[currentEntry] + i
      currentEntry = currentEntry - 1
    currentEntry = len(globalNumberOfOuterCells)-1
    for i in reversed(numberOfOuterCells[rank]):
      globalNumberOfOuterCells[currentEntry] = globalNumberOfOuterCells[currentEntry] + i
      currentEntry = currentEntry - 1

  pylab.clf()
  pylab.title( "Global cells" )
  x = pylab.arange(0, len(numberOfInnerLeafCells[0]), 1.0)
  pylab.plot(tTotal[0], globalNumberOfInnerLeafCells, 'o-',  markersize=10, color='#ff0000', label='#inner leaf cells' )
  pylab.plot(tTotal[0], globalNumberOfOuterLeafCells, '+-',  markersize=10, color='#00ff00', label='#outer leaf cells' )
  pylab.plot(tTotal[0], globalNumberOfInnerCells,     '.-',  markersize=10, color='#aa0066', label='#inner cells' )
  pylab.plot(tTotal[0], globalNumberOfOuterCells,     'x-',  markersize=10, color='#00aa66', label='#outer cells' )
  pylab.legend(fontsize=9, loc='upper left', framealpha=0.5)
  setGeneralPlotSettings()
  pylab.savefig( inputFileName + ".grid-overview.png" )
  pylab.savefig( inputFileName + ".grid-overview.pdf" )

  pylab.clf()
  pylab.title( "Local vertices" )
  startRank = 1
  if (numberOfRanks==1):
    startRank = 0
  for rank in range(startRank,numberOfRanks):
    x = pylab.arange(0, len(numberOfLocalVertices[rank]), 1.0)
    if len(tTotal[rank])==len(numberOfLocalVertices[rank]):
      pylab.plot(tTotal[rank], numberOfLocalVertices[rank], 'o',  color='#000000', alpha=AlphaValue, markersize=10)
  pylab.xlabel('t')
  pylab.savefig( inputFileName + ".local-vertices.png" )
  pylab.savefig( inputFileName + ".local-vertices.pdf" )


def plotWalltimeOverview(numberOfRanks,inputFileName): 
  pylab.clf()
  pylab.title( "Walltime" )
  pylab.ylabel( "time per grid sweep [t]=s" )
  x = pylab.arange(0, len(tTotal[0]), 1.0)
  pylab.plot(tTotal[0], tTraversal[0], '-',  markersize=10, color='#000066', label='time per traversal on global master' )
  startRank = 1
  if (numberOfRanks==1):
    startRank = 0
  for rank in range(startRank,numberOfRanks):
    x = pylab.arange(0, len(tTraversal[rank]), 1.0)
    if len(tTotal[rank])==len(tTraversal[rank]):
      pylab.plot(tTotal[rank], tTraversal[rank], 'o',  color='r', alpha=AlphaValue, markersize=10)
  setGeneralPlotSettings()
  pylab.savefig( inputFileName + ".walltime.png" )
  pylab.savefig( inputFileName + ".walltime.pdf" )





def plotStatisticsForRank(currentRank,numberOfRanks,inputFileName):
  pylab.clf()
  pylab.title( "Walltime" )
  pylab.plot(tTotal[0], tTraversal[0], '-',  markersize=10, color='#000066', label='time per traversal (global master)' )
  if len(tTotal[currentRank])==len(tTraversal[currentRank]):
    pylab.plot(tTotal[currentRank], tTraversal[currentRank], '-',  markersize=10, color='#550000', label='time per traversal (only inside local domain)' )
  else:
    print "WARNING: tTotal and tTraversal for rank " + str(currentRank) + " do not hold the same number of entries. Input file seems to be corrupted"
  startRank = 1
  if (numberOfRanks==1):
    startRank = 0
  for rank in range(startRank,numberOfRanks):
    if len(tTotal[rank])==len(tTraversal[rank]):
      pylab.plot(tTotal[rank], tTraversal[rank], 'o',  color='r', alpha=AlphaValue, markersize=10)
    else:
      print "WARNING: tTotal and tTraversal for rank " + str(rank) + " do not hold the same number of entries. Input file seems to be corrupted"
  setGeneralPlotSettings()
  pylab.savefig( inputFileName + ".walltime-rank-" + str(currentRank) + ".png" )
  pylab.savefig( inputFileName + ".walltime-rank-" + str(currentRank) + ".pdf" )

  pylab.clf()
  pylab.title( "Cells" )
  
  if len(tTotal[currentRank]) == len(numberOfLocalCells[currentRank]):
    pylab.plot(tTotal[currentRank], numberOfLocalCells[currentRank], '-',  markersize=10, color='b', label='local cells' )
  else:
    print "WARNING: arrays tTotal and numberOfLocalCells of rank " + str(currentRank) + " have different sizes. Input file might be corrupted"
  
  startRank = 1
  if (numberOfRanks==1):
    startRank = 0
  for rank in range(startRank,numberOfRanks):
    if len(tTotal[rank]) == len(numberOfLocalCells[rank]):
      pylab.plot(tTotal[rank], numberOfLocalCells[rank], 'o',  color='r', alpha=AlphaValue, markersize=10)
    else:
      print "WARNING: arrays tTotal and numberOfLocalCells of rank " + str(rank) + " have different sizes. Input file might be corrupted"
  setGeneralPlotSettings()
  pylab.savefig( inputFileName + ".local-cells-rank-" + str(currentRank) + ".png" )
  pylab.savefig( inputFileName + ".local-cells-rank-" + str(currentRank) + ".pdf" )
  
  totalTimeCalendar            = []
  joinTimeCalendar             = []
  boundaryTimeCalendar         = []
  synchronousHeapDataCalendar  = []
  asynchronousHeapDataCalendar = []
  centralElementCalendar       = []
  workerWaitTimeCalendar       = []

  totalTimeCPU                 = []
  joinTimeCPU                  = []
  boundaryTimeCPU              = []
  synchronousHeapDataCPU       = []
  asynchronousHeapDataCPU      = []
  centralElementCPU            = []

  ttotalTimeCalendar            = []
  tjoinTimeCalendar             = []
  tboundaryTimeCalendar         = []
  tsynchronousHeapDataCalendar  = []
  tasynchronousHeapDataCalendar = []
  tcentralElementCalendar       = []
  tworkerWaitTimeCalendar       = []

  ttotalTimeCPU                 = []
  tjoinTimeCPU                  = []
  tboundaryTimeCPU              = []
  tsynchronousHeapDataCPU       = []
  tasynchronousHeapDataCPU      = []
  tcentralElementCPU            = []

  try:
    inputFile = open( inputFileName,  "r" )
    for line in inputFile:
      if (re.search( "DefaultAnalyser::endIteration", line ) and re.search( "rank:" + str(currentRank) + " ", line ) and re.search( "t_total", line )):
        totalTimeCalendar.append( float( line.split( "(" )[-1].split(",")[0] ))
        totalTimeCPU.append(      float( line.split( "," )[-1].split(")")[0] ))
        
        ttotalTimeCalendar.append( float( line.strip().split( " " )[0] ))
        ttotalTimeCPU.append(      float( line.strip().split( " " )[0] ))
      if (re.search( "DefaultAnalyser::endReleaseOfJoinData", line ) and re.search( "rank:" + str(currentRank) + " ", line )):
        joinTimeCalendar.append( float( line.split( "=" )[1].split(",")[0] ))
        joinTimeCPU.append(      float( line.split( "=" )[2].strip() ))

        tjoinTimeCalendar.append( float( line.strip().split( " " )[0] ))
        tjoinTimeCPU.append(      float( line.strip().split( " " )[0] ))
      if (re.search( "DefaultAnalyser::endReleaseOfBoundaryData", line ) and re.search( "rank:" + str(currentRank) + " ", line )):
        boundaryTimeCalendar.append( float( line.split( "=" )[-2].split(",")[0] ))
        boundaryTimeCPU.append(      float( line.split( "=" )[-1].strip() ))

        tboundaryTimeCalendar.append( float( line.strip().split( " " )[0] ))
        tboundaryTimeCPU.append(      float( line.strip().split( " " )[0] ))
      if (re.search( "DefaultAnalyser::endToReleaseSynchronousHeapData", line ) and re.search( "rank:" + str(currentRank) + " ", line )):
        synchronousHeapDataCalendar.append( float( line.split( "=" )[-2].split(",")[0] ))      
        synchronousHeapDataCPU.append(      float( line.split( "=" )[-1].strip() ))

        tsynchronousHeapDataCalendar.append( float( line.strip().split( " " )[0] )    )  
        tsynchronousHeapDataCPU.append(      float( line.strip().split( " " )[0] ))
      if (re.search( "DefaultAnalyser::endToPrepareAsynchronousHeapDataExchange", line ) and re.search( "rank:" + str(currentRank) + " ", line )):
        asynchronousHeapDataCalendar.append( float( line.split( "=" )[-2].split(",")[0] ))      
        asynchronousHeapDataCPU.append(      float( line.split( "=" )[-1].strip() ))

        tasynchronousHeapDataCalendar.append( float( line.strip().split( " " )[0] ))     
        tasynchronousHeapDataCPU.append(      float( line.strip().split( " " )[0] ))
      if (re.search( "DefaultAnalyser::leaveCentralElementOfEnclosingSpacetree", line ) and re.search( "rank:" + str(currentRank) + " ", line )):
        centralElementCalendar.append( float( line.split( "(" )[-1].split(",")[0] ))
        centralElementCPU.append(      float( line.split( "," )[-1].split(")")[0] ))

        tcentralElementCalendar.append( float( line.strip().split( " " )[0] ))
        tcentralElementCPU.append(      float( line.strip().split( " " )[0] ))
      if (re.search( "DefaultAnalyser::dataWasNotReceivedFromWorker", line ) and re.search( "rank:" + str(currentRank) + " ", line )):
        workerWaitTimeCalendar.append( float( line.split( "for " )[2].split("s")[0] ))      
  
        tworkerWaitTimeCalendar.append( float( line.strip().split( " " )[0] ) )     
        
        
  except Exception as inst:
    print "failed to read " + inputFileName
    print inst

  pylab.clf()
  pylab.title( "Runtime profile (calendar time)" )

  if len(ttotalTimeCalendar)==len(totalTimeCalendar):
    pylab.plot(ttotalTimeCalendar, totalTimeCalendar, 'x-', label='total', color='#FF0000', markersize=10)
  else:
    print "ERROR: ttotalTimeCalendar has lenght " + str(len(ttotalTimeCalendar)) + " while totalTimeCalendar has lenght " + str(len(totalTimeCalendar))

  if len(tjoinTimeCalendar)==len(joinTimeCalendar):
    pylab.plot(tjoinTimeCalendar, joinTimeCalendar, 'o-', label='join', color='#00FF00', markersize=10)
  else:
    print "ERROR: tjoinTimeCalendar has lenght " + str(len(tjoinTimeCalendar)) + " while joinTimeCalendar has lenght " + str(len(joinTimeCalendar))

  if len(tboundaryTimeCalendar)==len(boundaryTimeCalendar):
    pylab.plot(tboundaryTimeCalendar, boundaryTimeCalendar, '^-', label='boundary exchange', color='#0000FF', markersize=10)
  else:
    print "ERROR: tboundaryTimeCalendar has lenght " + str(len(tboundaryTimeCalendar)) + " while boundaryTimeCalendar has lenght " + str(len(boundaryTimeCalendar))

  if len(tsynchronousHeapDataCalendar)==len(synchronousHeapDataCalendar):
    pylab.plot(tsynchronousHeapDataCalendar, synchronousHeapDataCalendar, '+-', label='synchronous heap', color='#aaaa00', markersize=10)
  else:
    print "ERROR: tsynchronousHeapDataCalendar has lenght " + str(len(tsynchronousHeapDataCalendar)) + " while synchronousHeapDataCalendar has lenght " + str(len(synchronousHeapDataCalendar))

  if len(tasynchronousHeapDataCalendar)==len(asynchronousHeapDataCalendar):
    pylab.plot(tasynchronousHeapDataCalendar, asynchronousHeapDataCalendar, '1-', label='asynchronous heap', color='#aa00aa', markersize=10)
  else:
    print "ERROR: tasynchronousHeapDataCalendar has lenght " + str(len(tasynchronousHeapDataCalendar)) + " while asynchronousHeapDataCalendar has lenght " + str(len(asynchronousHeapDataCalendar))

  #
  # If multiple ranks plot to one terminal, lines can be mixed up and parsing can fail. 
  #
  if (len(tcentralElementCalendar)==len(centralElementCalendar)):
    pylab.plot(tcentralElementCalendar, centralElementCalendar, 'v--', label='local elements traversal', color='#00aaaa', markersize=10)
  else:
    print "ERROR: tcentralElementCalendar has lenght " + str(len(tcentralElementCalendar)) + " while centralElementCalendar has lenght " + str(len(centralElementCalendar))

  pylab.plot(tworkerWaitTimeCalendar, workerWaitTimeCalendar, '3-', label='wait for worker', color='#667766', markersize=10)

  setGeneralPlotSettings()

  pylab.savefig( inputFileName + ".runtime-profile-calendar-rank-" + str(currentRank) + ".png" )
  pylab.savefig( inputFileName + ".runtime-profile-calendar-rank-" + str(currentRank) + ".pdf" )



  pylab.clf()
  pylab.title( "Runtime profile (cpu time)" )

  pylab.plot(ttotalTimeCPU, totalTimeCPU, 'x-', label='total', color='#FF0000', markersize=10)

  pylab.plot(tjoinTimeCPU, joinTimeCPU, 'o-', label='join', color='#00FF00', markersize=10)

  pylab.plot(tboundaryTimeCPU, boundaryTimeCPU, '^-', label='boundary exchange', color='#0000FF', markersize=10)

  pylab.plot(tsynchronousHeapDataCPU, synchronousHeapDataCPU, '+-', label='synchronous heap', color='#aaaa00', markersize=10)

  pylab.plot(tasynchronousHeapDataCPU, asynchronousHeapDataCPU, '1-', label='asynchronous heap', color='#aa00aa', markersize=10)

  if (len(tcentralElementCPU)==len(centralElementCPU)):
    pylab.plot(tcentralElementCPU, centralElementCPU, 'v--', label='local elements traversal', color='#00aaaa', markersize=10)
  else:
    print "- error in centralElementCPU"

  setGeneralPlotSettings()
  pylab.savefig( inputFileName + ".runtime-profile-cpu-rank-" + str(currentRank) + ".png" )
  pylab.savefig( inputFileName + ".runtime-profile-cpu-rank-" + str(currentRank) + ".pdf" )

  
def computeVolumesOverlapsWork(numberOfRanks,volume,offset,dim,domainoffset,domainsize,parents):
  print "compute volumes, overlaps and work ",
  volumes = [0.0 for x in offset]
  for i in range(0,numberOfRanks):
    print ".",
    if dim==2:
      volumes[i] = volume[i][0]*volume[i][1]
    else:
      volumes[i] = volume[i][0]*volume[i][1]*volume[i][2]

  overlaps = [0.0 for x in offset]
  for i in range(0,numberOfRanks):
    print ".",
    overlaps[i] = 1.0
    for d in range(0,dim):
      left  = max(offset[i][d],float(domainoffset[d]))
      right = min(offset[i][d]+volume[i][d],float(domainoffset[d])+float(domainsize[d]))
      delta = right-left
      if delta<0:
        print "ERROR for rank " + str(i) +  \
              ": region=" + str(offset[i]) + "x" + str(volume[i]) + \
              ", bonding box=" + str(domainoffset) + "x" + str(domainsize) + \
              ", delta= " + str(delta) 
      overlaps[i] = overlaps[i] * delta
    if overlaps[i]>volume[i]:
      print "ERROR for rank " + str(i) +  \
            ": region=" + str(offset[i]) + "x" + str(volume[i]) + \
            ", bounding box=" + str(domainoffset) + "x" + str(domainsize) + \
            ", volume= " + str(volumes[i]) + ", overlaps=" + str(overlaps[i]) 
  work = [x for x in overlaps]
  for i in range(1,numberOfRanks):
    if work[ parents[i] ]<overlaps[i]:
      print "WARNING: work of rank " + str(parents[i]) + " will become negative as overlap of rank " + str(i) + " equals " + str(overlaps[i]) 
    work[ parents[i] ] = work[ parents[i] ] - overlaps[i]
  work[0] = 0.0     # can become negative

  print " done "
  return (volumes,overlaps,work)
  
  
  
  
def plotWorkloadAndResponsibilityDistribution(numberOfRanks,volumes,overlaps,work,outputFileName):
 print "plot workload distribution ",
 pylab.clf()
 ranks = [x for x in range(0,numberOfRanks)]

 print ".",
 _markevery = numberOfRanks
 if numberOfRanks>11:
  _markevery = _markevery / 11
 pylab.plot(ranks, volumes, '-o', label="area of responsibility", markevery=_markevery, color='#000000', markersize=10)

 print ".",
 _markevery = numberOfRanks
 if numberOfRanks>11:
   _markevery = _markevery / 11
 pylab.plot(ranks, overlaps, '-s', label="area of responsibility $\cap \Omega$", markevery=_markevery, color='#ff0000', markersize=10)

 print ".",
 pylab.fill_between(ranks, work, color='#0000bb', alpha=0.4)

 maxLocalWorkToGetExtraLabel = 0
 ranksWithLabels = len(work)
 while ranksWithLabels > numberOfRanks/10:
   maxLocalWorkToGetExtraLabel = 0.5 * maxLocalWorkToGetExtraLabel + 0.5 * max(work)
   ranksWithLabels = sum( i>maxLocalWorkToGetExtraLabel for i in work )
 
 for i in range(0,numberOfRanks):
   if work[i]>maxLocalWorkToGetExtraLabel:
     pylab.text(i,volumes[i]+10,str(i))
     pylab.plot([i,i], [0, volumes[i]], '--', color="#000000") 

 ranksWithZeroResponsibility = 0
 while volumes[ranksWithZeroResponsibility]==0.0:
   ranksWithZeroResponsibility = ranksWithZeroResponsibility + 1
 if (ranksWithZeroResponsibility>2):
   pylab.text(ranksWithZeroResponsibility,volumes[ranksWithZeroResponsibility]+20,"ranks_per_node")
   pylab.plot([ranksWithZeroResponsibility,ranksWithZeroResponsibility], [0, volumes[ranksWithZeroResponsibility]], '-', color="#000000") 
 
 try:
   pylab.legend(fontsize=9, framealpha=0.5)
 except:
   # old pylab version
   l = pylab.legend(prop={'size':9})
   l.get_frame().set_alpha(0.5)
 ax = pylab.gca()
 ax.autoscale_view()
 ax.set_yscale('symlog', basey=10)

 #if numberOfRanks>32:
 # pylab.figure(figsize=(numberOfRanks/10,4))

 setGeneralPlotSettings()
 pylab.xlabel('rank')  
 pylab.ylabel('$\Omega $')
 pylab.savefig( outputFileName + "-symlog.work-distribution.png" )
 pylab.savefig( outputFileName + "-symlog.work-distribution.pdf" )
 switchToLargePlot()
 pylab.savefig( outputFileName + "-symlog.work-distribution.large.png" )
 pylab.savefig( outputFileName + "-symlog.work-distribution.large.pdf" )
 switchBackToStandardPlot()

 ax.set_yscale('symlog', basey=10)

 #if numberOfRanks>32:
 # pylab.figure(figsize=(numberOfRanks/10,4))

 ax.set_yscale('linear' )
 setGeneralPlotSettings()
 pylab.xlabel('rank')  
 pylab.ylabel('$\Omega $')
 pylab.savefig( outputFileName + ".work-distribution.png" )
 pylab.savefig( outputFileName + ".work-distribution.pdf" )
 switchToLargePlot()
 pylab.savefig( outputFileName + ".work-distribution.large.png" )
 pylab.savefig( outputFileName + ".work-distribution.large.pdf" )
 switchBackToStandardPlot()
 print "done"

 
 
def plot2dDomainDecompositionOnLevel(l,numberOfRanks,domainoffset,domainsize,offset,volume,levels,outputFileName):
  Colors = [ "#ddeedd", "#ccbbcc", "#aabbaa", "#887788" ]
  print "plot domain decomposition on level " + str(l),
  pylab.clf()
  pylab.figure(figsize=(float(domainsize[0]),float(domainsize[1])))
  colorCounter = 0
  for i in range(0,numberOfRanks):
    if levels[i]==l:
      print ".",
      
      pylab.gca().add_patch(pylab.Rectangle(offset[i]+(volume[i][0]*0.01,volume[i][0]*0.01), volume[i][0]*0.98, volume[i][1]*0.98, color=Colors[colorCounter % len(Colors)]))
      
      chosenFontSize = int(round(float(domainsize[0])*12.0/(l)))
      colorCounter   = colorCounter + 1
      
      pylab.text(
        offset[i][0] + volume[i][0]/2,
        offset[i][1] + volume[i][1]/2,
        str(i), fontsize=chosenFontSize
      )
  print ".",
  pylab.xlim( float(domainoffset[0]), float(domainoffset[0])+float(domainsize[0]) )
  pylab.ylim( float(domainoffset[1]), float(domainoffset[1])+float(domainsize[1]) )
  pylab.savefig( outputFileName + ".level" + str(l) + ".png" )
  pylab.savefig( outputFileName + ".level" + str(l) + ".pdf" )
  print "done"
