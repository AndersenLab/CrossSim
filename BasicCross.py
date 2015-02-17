#!/usr/bin/env python
# encoding: utf-8
"""
WormIndividual.py

Zifan Xiang
Copyright (c) 2015 Northwestern University. All rights reserved.
"""

from Chromosomes import *
from CrossUtils import *
from Individual import *
from WormUtils import *
import numpy as np
import operator
import os.path
import random
import sys

# Crosses each worm in set A with worm B until the limit number of offspring that keep the parent segment at the desired location has been met 
def backCross(diploidASet, diploidB, physLoc, chromNumber, allele):
  generation = []
  loc = Chromosome.getLoc(physLoc, chromNumber)
  
  for diploidA in diploidASet:  
    curDiploid = diploidA.mate(diploidB)[0]
    
    while (curDiploid.chromosome_set[0][chromNumber].getParentAtLocation(loc) != allele and curDiploid.chromosome_set[1][chromNumber].getParentAtLocation(loc) != allele):                
      curDiploid = diploidA.mate(diploidB)[0]
      
    generation.append(curDiploid)
    
  return generation

def selfCross(diploidSet):
  generation = []
  loc = Chromosome.getLoc(physLoc, chromNumber)

  for diploid in diploidSet:
    curChild = diploid.mate(diploid)[0]
    generation.append(curChild)

  return generation

def randomCross(diploidSet):
  generation = []
  length = len(diploidSet)
  random.seed()

  for i in range(length):
    firstRandIndex = random.randint(0, length - 1)
    secondRandIndex = random.randint(0, length - 1)

    while (secondRandIndex == firstRandIndex):
      secondRandIndex = random.randint(0, length - 1)

    curChild = diploidSet[firstRandIndex].mate(diploidSet[secondRandIndex])[0]
    generation.append(curChild)
 
  return generation

def circularCross(diploidSet):
  generation = []
  length = len(diploidSet)

  for i in range(length):
    curChild = diploidSet[i].mate(diploidSet[(i + 1) % length])[0]
    generation.append(curChild)

  return generation

def randomCrossEqualContribution(diploidSet):
  generation = []
  length = len(diploidSet)
  random.seed()
  firstParentIndices = range(length)
  secondParentIndices = range(length)
  random.shuffle(firstParentIndices)
  random.shuffle(secondParentIndices)

  for i in range(len(firstParentIndices)):
    if (firstParentIndices[i] == secondParentIndices[i]):
      j = i + 1
      temp = secondParentIndices[j % len(secondParentIndices)]
      secondParentIndices[j % len(secondParentIndices)] = secondParentIndices[i]
      secondParentIndices[i] = temp

  for i in range(len(firstParentIndices)):
    curChild = diploidSet[firstParentIndices[i]].mate(diploidSet[secondParentIndices[i]])[0]
    generation.append(curChild)

  return generation

def circularPairCross(diploidSet):
  generation = []
  length = len(diploidSet)
  
  for i in range(0, length, 2):
    for j in range(2):
      curChild = diploidSet[i].mate(diploidSet[(i + 1) % length])[0]
      generation.append(curChild)
    
    # Shuffling the order of the next generation
    generation = generation[len(generation) - 1:] + generation[1:]

  return generation

def randomPairCross(diploidSet):
  generation = []
  length = len(diploidSet)
  random.seed()
  parentIndices = range(length)
  random.shuffle(parentIndices)

  for i in range(length):
    curIndex = random.randrange(0, length - 2, 2)
    curChild = diploidSet[parentIndices[curIndex]].mate(diploidSet[parentIndices[curIndex + 1]])[0]
    generation.append(curChild)

  return generation

def randomPairCrossEqualContribution(diploidSet):
  generation = []
  length = len(diploidSet)
  parentIndices = range(length)
  random.shuffle(parentIndices)

  for i in range(length / 2):
    for j in range(2):
      curChild = diploidSet[parentIndices[i]].mate(diploidSet[parentIndices[i + 1]])[0]
      generation.append(curChild)

  return generation

def inbreedingAvoidanceCross(diploidSet):
  length = len(diploidSet)
  generation = []

  for k in range(2):
    for i in range(length / 2):
      curChild = diploidSet[i].mate(diploidSet[i + 1])[0]
      generation.append(curChild)

  return generation

def writeGeneralStatistics(crossNumber, physLoc, diploidSet, targetChrom, targetName, statFile):
  indNumber = 1;
  genLoc = Chromosome.getLoc(physLoc, chromNumber)
  
  for diploid in diploidSet:
    totalSelected = 0
    curIntervals = []
  
    for chrSet in diploid.chromosome_set:
      percent = chrSet[targetChrom].getPercentageOfParent(targetName)
      totalSelected += percent
        
      if chrSet[chromNumber].getParentAtLocation(genLoc) == targetName:
        curIntervals.append(chrSet[chromNumber].physicalLocsOfInterval(genLoc, chromNumber))
            
    totalLower = 0;
    totalUpper = 0;
    for interval in curIntervals:
      totalLower += interval[0];
      totalUpper += interval[1]
      
    perGenome = diploid.getPercentageOfGenome(targetName, chromNumber) # Percent of genome not including the selected chromosome
    avgLower = totalLower / len(curIntervals)
    avgUpper = totalUpper / len(curIntervals)
    avgSelected = totalSelected / 2
    statFile.write('%d,%d,%d,%d,%f,%f,%d,%d\n' % (crossNumber, indNumber, targetChrom + 1, physLoc, avgSelected, perGenome, avgLower, avgUpper))
    indNumber += 1

#Takes in the general statistics file and writes summary statistics on a particular generation of the back cross    
def writeGroupSegments(fileName, diploidSet):
  f = open(fileName, 'wb')
  f.write('Individual,Set,Chromosome,Type,Left Position,Right Position\n')
  t = 1 #Counter for the different individuals that have to printed within each cross 
  for diploid in diploidSet:
    l = 1 #Counter for the two sets of chromosomes
    for chrSet in diploid.chromosome_set:
      j = 1
      for chromosome in chrSet:
        for i in range(len(chromosome.segments)):
          parent = chromosome.segments[i][1]
          leftLoc = Chromosome.getPhysDistanceFromLoc(chromosome.segments[i][0], j - 1) 
          rightLoc = Chromosome.getPhysDistanceFromLoc(1, j - 1) ;
          
          if (i + 1 != len(chromosome.segments)):
            rightLoc = Chromosome.getPhysDistanceFromLoc(chromosome.segments[i + 1][0], j - 1) 
        
          f.write('%d,%d,%d,%s,%d,%d\n' % (t, l, j, parent, leftLoc, rightLoc))
        j += 1
      l += 1
    t += 1
    
  f.close()

#num represents if the desired intervals are the lower or higher intervals, num = 0/1
def separatePhysicalInterval(selectedPhysInterval, num):
  toReturn = [0 for x in range(len(selectedPhysInterval))]
  
  for i in range(len(selectedPhysInterval)):
    toReturn[i] = selectedPhysInterval[i][num]
    
  return toReturn

def bucketPhysicalIntervals(selectedPhysInterval, low, high, bucketSize):
  size = int(high - low) / bucketSize
  buckets = [0 for i in range(size)]
  
  if (len(buckets) == 0):
    buckets = [0]
  
  for i in range(len(selectedPhysInterval)):
    bucketPos = int((selectedPhysInterval[i] - low) / bucketSize)
    
    if bucketPos == len(buckets):
      bucketPos -= 1

    buckets[bucketPos] += 1

  return buckets

def putIntervalsIntoBuckets(numCross, chromNumber, physLoc, physInterval, bucketSize, numSampled, filePath):    
  filePath.write('%d,%d,%d,%d,%d,' % (numCross, chromNumber + 1, physLoc, bucketSize, numSampled))
  
  for i in range(2):
    sepPhysicalInterval = separatePhysicalInterval(physInterval, i)
    low = min(sepPhysicalInterval)
    high = max(sepPhysicalInterval)
    toWrite = bucketPhysicalIntervals(sepPhysicalInterval, low, high, bucketSize)
    numUnique = 0;
    
    for bucket in toWrite:
      if bucket != 0:
        numUnique += 1;
        
    filePath.write('%d,%d,%d' % (low, high, numUnique))
    
    if (i == 0): 
      filePath.write(',')
      
  filePath.write('\n')

#Selects a random subset from a set 
def selectRandomSubset(diploidSet, numSelect):
  length = len(diploidSet)
  randomIndices = [];
  toReturnSet = []
  random.seed()
  
  if numSelect > length:
    raise ValueError, "The number to be selected cannot be greater than the total set, repeats would be produced"
    
  while (len(randomIndices) < numSelect):
    curIndex = random.randint(0, length - 1)
    
    if (not(curIndex in randomIndices)):
      randomIndices.append(curIndex)
  
  for i in randomIndices:
    toReturnSet.append(diploidSet[i])
  
  return toReturnSet

def writeAverageBinSize(crossNumber, diploidSet, writeFile):
  averageBinSizes = []

  for diploid in diploidSet:
    averageBinSizes.append(diploid.getAverageBinGeneticSize())

  statisticsWrapper = np.array(averageBinSizes)
  writeFile.write('%d,%f,%f\n' % (crossNumber, np.mean(statisticsWrapper), np.std(statisticsWrapper)))

#Simulates a back crosses in which a particular base pair allele is held 
def backCrossSimulation(physLoc, chromNumber, crossNumber, numIndividuals, bucketSize, numRandomSelect, numIter, crossOption):
  #Opens a file that contains info about the general statistics (percentage of genome, percentage of selected chromosome) of the cross simulation
  if os.path.isfile('general_statistics_%d_%d.csv' % (physLoc, chromNumber + 1)):
    g = open('general_statistics_%d_%d.csv' % (physLoc, chromNumber + 1), 'a')
  else:
    g = open('general_statistics_%d_%d.csv' % (physLoc, chromNumber + 1), 'wb')
    g.write('Number of Back Crosses,Individual Number,Selected Chromosome,Selected Base Pair,Percent Selected Chromosome,Percent Genome,Left Physical Loc,Right Physical Loc\n')
  
  #Opens a file that contains the info about the number of unique intervals  
  if os.path.isfile('buckets_%d_%d_%d_%d.csv' % (physLoc, chromNumber + 1, bucketSize, numRandomSelect)):
    h = open('buckets_%d_%d_%d_%d.csv' % (physLoc, chromNumber + 1, bucketSize, numRandomSelect), 'a')
  else:
    h = open('buckets_%d_%d_%d_%d.csv' % (physLoc, chromNumber + 1, bucketSize, numRandomSelect), 'wb')
    h.write('Number of Back Crosses,Selected Chromosome,Selected Base Pair,Bucket Size,Number Sampled,Minimum Left Base Pair, Maximum Left Base Pair,Number Left Unique Buckets,Minimum Right Base Pair, Maximum Right Base Pair,Number Right Unique Buckets\n')
  
  if os.path.isfile('average_genetic_bin_size_%d_%d.csv' % (physLoc, chromNumber + 1)):
    binSizeFile = open('average-genetic_bin_size_%d_%d.csv' % (physLoc, chromNumber + 1), 'a')
  else:
    binSizeFile = open('average_genetic_bin_size_%d_%d.csv' % (physLoc, chromNumber + 1), 'wb')   
    binSizeFile.write('Number of Back Crosses,Average Bin Size,Bin Size Standard Deviation\n')

  #Runs through the number of crosses specified and makes the individuals
  diploidSet = generateHetero(numIndividuals)
        
  Bparent = Diploid(name = "B", newChr = 6)
  targetNameDip = "A"
  genLoc = Chromosome.getLoc(physLoc, chromNumber)

  for k in range(crossNumber):
    #TODO(zifanxiang): This if else ladder is disgusting
    if crossOption == 1:
      diploidSet = backCross(diploidSet, Bparent, physLoc, chromNumber, targetNameDip)
    elif crossOption == 2:
      diploidSet = selfCross(diploidSet)
    elif crossOption == 3:
      diploidSet = circularCross(diploidSet)
    elif crossOption == 4:
      diploidSet = circularPairCross(diploidSet)
    elif crossOption == 5:
      diploidSet = inbreedingAvoidanceCross(diploidSet)
    elif crossOption == 6:
      diploidSet = randomCross(diploidSet)
    elif crossOption == 7:
      diploidSet = randomCrossEqualContribution(diploidSet)
    elif crossOption == 8:
      diploidSet = randomPairCross(diploidSet)
    elif crossOption == 9:
      diploidSet = randomPairCrossEqualContribution(diploidSet)
    else:
      diploidSet = backCross(diploidSet, Bparent, physLoc, chromNumber, targetNameDip)
    
    #writeGeneralStatistics(k + 1, physLoc, diploidSet, chromNumber, targetNameDip, g)
    writeAverageBinSize(k + 1, diploidSet, binSizeFile)
    
    for i in range(numIter):
      physIntervals = []
      sampleSet = selectRandomSubset(diploidSet, numRandomSelect);
      
      for diploid in sampleSet:
        for chrSet in diploid.chromosome_set:
          if chrSet[chromNumber].getParentAtLocation(genLoc) == targetNameDip:
            physIntervals.append(chrSet[chromNumber].physicalLocsOfInterval(genLoc, chromNumber))
    
      #putIntervalsIntoBuckets(k + 1, chromNumber, physLoc, physIntervals, bucketSize, numRandomSelect, h)
      
    # Format of the output files is as follows: Number of Crosses_ Number Of Individuals per Cross _ Target Chromosome _ Physical Location on the Target Chromosome
    fileName = "%d_%d_%d_%d_crossConfig.csv" % (k + 1, numIndividuals, chromNumber + 1, physLoc)
  
    truncAparentSet = selectRandomSubset(diploidSet, numRandomSelect)    
    writeGroupSegments(fileName, truncAparentSet)

  g.close()

#Parameters: physLoc chromNumber numCrosses numIndividuals bucketSize numRandomSelect numIter
if __name__ == '__main__':
  physLoc = int(sys.argv[1]) #Physical location that must be held
  chromNumber = int(sys.argv[2]) - 1;
  numCrosses = int(sys.argv[3])
  numIndividuals = int(sys.argv[4])
  bucketSize = int(sys.argv[5])
  numRandomSelect = int(sys.argv[6])
  numIter = int(sys.argv[7])
  crossOption = int(sys.argv[8])
  
  for i in range(numIter):
    backCrossSimulation(physLoc, chromNumber, numCrosses, numIndividuals, bucketSize, numRandomSelect, numIter, crossOption)