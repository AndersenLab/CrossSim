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

from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

import operator
import os.path
import random
import sys
import threading

crossDesign = ["self", "circular", "circular_pair", "interbreed_avoidance", "random", "random_equal", "random_pair", "random_pair_equal"]

def selfCross(diploidSet):
  generation = []

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
    generation = generation[len(generation) - 1:] + generation[0:len(generation) - 1];

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

# Uses Individual.py's getExpectedBinGeneticSizes method, which calculates the average bin sizes (interval) in
# centimorgans (converts the genetic location between 0 to 1 through the centimorgan lengths of the chromosomes)
def calcExpectedBinSize(diploidSet):
  averageBinSizes = []

  for diploid in diploidSet:
    averageGeneticBinSizes = diploid.getExpectedBinGeneticSizes()

    for size in averageGeneticBinSizes:
      averageBinSizes.append(size)

  return averageBinSizes

# Calculates the percent deviation from 50% of a chromosome at a particular location
def calcGeneticDrift(diploidSet, chromNumber, randomMarkerLocs, targetAllele):
  geneticDrifts = []

  for randomMarkerLoc in randomMarkerLocs:
    curDrifts = []
    loc = Chromosome.getLoc(randomMarkerLoc, chromNumber)
    for diploid in diploidSet:
      count = diploid.getNumAllele(chromNumber, targetAllele, loc)
      curDrifts.append(float(count) / 2)
        
    geneticDrifts.append(abs(np.mean(np.array(curDrifts)) - 0.5))
    
  return geneticDrifts

# Calculates the number breakpoints present in all of the chromosomes in the generation
def calcActualBreakpoints(diploidSet):
  numBreakpoints = 0

  for diploid in diploidSet:
    numBreakpoints = numBreakpoints + diploid.getNumBreakpoints()

  return numBreakpoints

def generateNextGeneration(crossOption, diploidSet):
  if crossOption == 1:
    diploidSet = selfCross(diploidSet)
  elif crossOption == 2:
    diploidSet = circularCross(diploidSet)
  elif crossOption == 3:
    diploidSet = circularPairCross(diploidSet)
  elif crossOption == 4:
    diploidSet = inbreedingAvoidanceCross(diploidSet)
  elif crossOption == 5:
    diploidSet = randomCross(diploidSet)
  elif crossOption == 6:
    diploidSet = randomCrossEqualContribution(diploidSet)
  elif crossOption == 7:
    diploidSet = randomPairCross(diploidSet)
  elif crossOption == 8:
    diploidSet = randomPairCrossEqualContribution(diploidSet)

  return diploidSet

def simulateRecombinantInbredCross(crossBinSizes, geneticDrifts, mapExpansions, numSelfing, k, crossOption, chromNumber, randomMarkerLocs, targetAllele, z):
  diploidSet = generateHeterozygotes(numIndividuals)
  diploidSet = generateNextGeneration(crossOption, diploidSet)
  print diploidSet[0].chromosome_set[0][0].interference

  for j in range(k):
    diploidSet = generateNextGeneration(crossOption, diploidSet)

    for j in range(numSelfing):
      diploidSet = selfCross(diploidSet)

      curAvgBinSizes = np.array(calcExpectedBinSize(diploidSet))
      curAvgBinSize = np.mean(curAvgBinSizes)
      curAvgBinSizeStd = np.std(curAvgBinSizes)
      curGeneticDrifts = np.array(calcGeneticDrift(diploidSet, chromNumber, randomMarkerLocs, targetAllele))
      curGeneticDrift = np.mean(curGeneticDrifts)
      curGeneticDriftStd = np.std(curGeneticDrifts)
      expectedNumBreakpoints = (float(numIndividuals) / 4) * ((k / 2) + (pow(2, numSelfing) - 1) / (pow(2, numSelfing))) 
      crossBinSizes[z].append(curAvgBinSize)
      geneticDrifts[z].append(np.mean(np.array(calcGeneticDrift(diploidSet, chromNumber, randomMarkerLocs, targetAllele))))
      mapExpansions[z].append(calcActualBreakpoints(diploidSet) / expectedNumBreakpoints)

def recombinantInbredLinesSimulation(numIndividuals, numCrosses, numSelfing, numIter, crossOption):
  if os.path.isfile('expected_genetic_bin_size_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1])):
    binSizeFile = open('expected_genetic_bin_size_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1]), 'a')
  else:
    binSizeFile = open('expected_genetic_bin_size_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1]), 'wb')   
    binSizeFile.write('Number of Back Crosses,Average Bin Size,Bin Size Standard Deviation\n')

  if os.path.isfile('genetic_drift_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1])):
    geneticDriftFile = open('genetic_drift_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1]))
  else:
    geneticDriftFile = open('genetic_drift_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1]), 'wb')
    geneticDriftFile.write('Number of Back Crosses,Genetic Drift,Standard Deviation\n')

  if os.path.isfile('map_expansion_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1])):
    mapExpansionFile = open('map_expansion_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1]), 'a')
  else:
    mapExpansionFile = open('map_expansion_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1]), 'wb')   
    mapExpansionFile.write('Number of Back Crosses,Map Expansion Ratio,Standard Deviation\n')

  start = datetime.now()
  targetAllele = 'A' if random.randint(0, 1)  == 0 else 'B'
  crossBinSizes = [[] for x in range(4, numCrosses + 1, 2)]
  geneticDrifts = [[] for x in range(4, numCrosses + 1, 2)]
  mapExpansions = [[] for x in range(4, numCrosses + 1, 2)]
  randomMarkerLocs = []
  chromNumber = random.randint(0, 5)

  for i in range(1):
    randomMarkerLocs.append(random.randint(0, chromosome_phys_max[chromNumber] - 1))

  z = 0
  # Outer loop runs through the different number of crosses, k represents the number of crosses in the RIAIL
  #for k in range(4, numCrosses + 1, 2):
    #threads = []
    #for i in range(numIter):
    #  t = threading.Thread(target=simulateRecombinantInbredCross, args=(crossBinSizes, geneticDrifts, mapExpansions, numSelfing, k, crossOption, chromNumber, randomMarkerLocs, targetAllele, z))
    #  threads.append(t)
    #  t.start();
    #  print "Thread number %d with num crosses %d has started" % (i, k)

    #for thread in threads:
    #  thread.join()

  #  z = z + 1

  for k in range(4, numCrosses + 1, 2):
    print k
    print '\n'
    for i in range(numIter):
     simulateRecombinantInbredCross(crossBinSizes, geneticDrifts, mapExpansions, numSelfing, k, crossOption, chromNumber, randomMarkerLocs, targetAllele, z)

    binSizeFile.write('%d,%f,%f\n' % (k, np.mean(np.array(crossBinSizes[z])), np.std(np.array(crossBinSizes[z]))))
    geneticDriftFile.write('%d,%f,%f\n' % (k, np.mean(np.array(geneticDrifts[z])), np.std(np.array(geneticDrifts[z]))))
    mapExpansionFile.write('%d,%f,%f\n' % (k, np.mean(np.array(mapExpansions[z])), np.std(np.array(mapExpansions[z]))))
    print i
    print '\n'

    z = z + 1

  end = datetime.now()
  diff = end - start;
  print diff.seconds

  plt.boxplot(geneticDrifts)
  plt.xlabel("Number of Crosses")
  plt.ylabel("Deviation from 50 Percent")
  plt.title("Genetic Drift")
  plt.xticks(np.arange(1, len(range(4, numCrosses + 1))), range(4, numCrosses + 1, 2))
  plt.savefig("geneticDrift_%d_%s_%d_%d_%d.png" % (numIndividuals, crossDesign[crossOption - 1], numIter, numCrosses, numSelfing), bbox_inches='tight')
  plt.show()

  plt.boxplot(crossBinSizes)
  plt.xlabel("Number of Crosses")
  plt.ylabel("Expected Bin Size (cM)")
  plt.title("Expected Bin Size")
  plt.xticks(np.arange(1, len(range(4, numCrosses + 1))), range(4, numCrosses + 1, 2))
  plt.savefig("expectedBinSize_%d_%s_%d_%d_%d.png" % (numIndividuals, crossDesign[crossOption - 1], numIter, numCrosses, numSelfing), bbox_inches='tight')
  plt.show()

  plt.boxplot(mapExpansions)
  plt.xlabel("Number of Crosses")
  plt.ylabel("Map Expansion Relative to Infinite Population Case")
  plt.title("Map Expansion")
  plt.xticks(np.arange(1, len(range(4, numCrosses + 1))), range(4, numCrosses + 1, 2))
  plt.savefig("mapExpansion_%d_%s_%d_%d_%d.png" % (numIndividuals, crossDesign[crossOption - 1], numIter, numCrosses, numSelfing), bbox_inches='tight')
  plt.show()

# Parameters for RecombinantInbredCross.py: 
# numIndividuals - the number of individuals per generation
# numCrosses - the max number of inbreeding generations to complete (goes in multiples of 2 from 4 up to the max)
# numSelfing - the number of selfing generations after the inbreeding generations
# numIter - the number of iterations to run each RIAIL
# crossOption - the type of crossing done during the inbreeding phase
if __name__ == '__main__':
  numIndividuals = int(sys.argv[1])
  numCrosses = int(sys.argv[2])
  numSelfing = int(sys.argv[3])
  numIter = int(sys.argv[4])
  crossOption = int(sys.argv[5])

  recombinantInbredLinesSimulation(numIndividuals, numCrosses, numSelfing, numIter, crossOption)