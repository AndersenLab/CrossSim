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
import matplotlib.pyplot as plt

import operator
import os.path
import random
import sys

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

# Uses Individual.py's getAverageBinGeneticSizes method, which calculates the average bin sizes (interval) in
# centimorgans (converts the genetic location between 0 to 1 through the centimorgan lengths of the chromosomes)
def calcExpectedBinSize(diploidSet):
  averageBinSizes = []

  for diploid in diploidSet:
    averageGeneticBinSizes = diploid.getExpectedBinGeneticSizes()

    for size in averageGeneticBinSizes:
      averageBinSizes.append(size)

  return averageBinSizes

def calcGeneticDrift(chromNumber, randomMarkerLocs, diploidSet, targetAllele):
  geneticDrifts = []

  for randomMarkerLoc in randomMarkerLocs:
    curDrifts = []
    loc = Chromosome.getLoc(randomMarkerLoc, chromNumber)
    for diploid in diploidSet:
      count = 0
    
      for chrSet in diploid.chromosome_set:
        if chrSet[chromNumber].getParentAtLocation(loc) == targetAllele:
          count += 1

      curDrifts.append(float(count) / 2)
        
    geneticDrifts.append(abs(np.mean(np.array(curDrifts)) - 0.5))
    
  return geneticDrifts

def recombinantInbredLinesSimulation(numIndividuals, numCrosses, numSelfing, numIter, crossOption):
  if os.path.isfile('average_genetic_bin_size_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1])):
    binSizeFile = open('average_genetic_bin_size_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1]), 'a')
  else:
    binSizeFile = open('average_genetic_bin_size_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1]), 'wb')   
    binSizeFile.write('Number of Back Crosses,Average Bin Size,Bin Size Standard Deviation\n')

  if os.path.isfile('genetic_drift_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1])):
    geneticDriftFile = open('genetic_drift_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1]))
  else:
    geneticDriftFile = open('genetic_drift_%d_%d_%d_%d_%s.csv' % (numIndividuals, numCrosses, numSelfing, numIter, crossDesign[crossOption - 1]), 'wb')
    geneticDriftFile.write('Number of Back Crosses,A Genetic Drift, B Genetic Drift,A Stand Dev,B Stand Dev\n')

  targetAllele = 'A' if random.randint(0, 1)  == 0 else 'B'
  crossBinSizes = [[] for x in range(4, numCrosses + 1, 2)]
  geneticDrifts = [[] for x in range(4, numCrosses + 1, 2)]
  randomMarkerLocs = []
  chromNumber = random.randint(0, 5)

  for i in range(10):
    randomMarkerLocs.append(random.randint(0, chromosome_phys_max[chromNumber] - 1))

  z = 0
  for k in range(4, numCrosses + 1, 2):
    for i in range(numIter):
      diploidSet = generateHeterozygotes(numIndividuals)
  
      for j in range(k):
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

      for j in range(numSelfing):
        diploidSet = selfCross(diploidSet)

      curAvgBinSizes = calcExpectedBinSize(diploidSet)
      crossBinSizes[z].append(np.mean(np.array(curAvgBinSizes)))
      geneticDrifts[z].append(np.mean(np.array(calcGeneticDrift(chromNumber, randomMarkerLocs, diploidSet, targetAllele))))

    z = z + 1

  plt.boxplot(geneticDrifts)
  plt.xlabel("Number of Crosses")
  plt.ylabel("Deviation from 50 Percent")
  plt.xticks(range(4, numCrosses + 1, 2))
  plt.show()

  plt.boxplot(crossBinSizes)
  plt.xlabel("Number of Crosses")
  plt.ylabel("Bin Size (cM)")
  plt.show()

if __name__ == '__main__':
  numIndividuals = int(sys.argv[1])
  numCrosses = int(sys.argv[2])
  numSelfing = int(sys.argv[3])
  numIter = int(sys.argv[4])
  crossOption = int(sys.argv[5])

  recombinantInbredLinesSimulation(numIndividuals, numCrosses, numSelfing, numIter, crossOption)