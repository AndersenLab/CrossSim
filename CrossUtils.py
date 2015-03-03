#!/usr/bin/env python
# encoding: utf-8
"""
CrossUtils.py
"""

from Chromosomes import *
from Individual import*

#Gets the average percentages of the target chromsome and genome from a set of individuals  
def averagePercentages(diploidSet, targetChrom, targetName):
  length = len(diploidSet)
  totalSelected = 0
  totalGenome = 0
    
  for diploid in diploidSet:
    for chrSet in diploid.chromosome_set:
      totalSelected += chrSet[targetChrom].getPercentageOfParent(targetName)
        
    totalGenome += diploid.getPercentageOfGenome(targetName)
    
  return [totalSelected / (length * 2), totalGenome / length]

#Calculates the average left and right intervals from a set of physical intervals       
def calculateAveragePhysicalIntervals(physIntervals, physLoc, chromNumber, filePath):
  lowerIntervalSums = 0
  upperIntervalSums = 0
  lowerIntervalAverages = 0
  upperIntervalAverages = 0
    
  for x in range(len(physIntervals)):
    interval = physIntervals[x]
    lowerIntervalSums += interval[0]
    upperIntervalSums += interval[1]

    lowerIntervalAverages = lowerIntervalSums / len(physIntervals)
    upperIntervalAverages = upperIntervalSums / len(physIntervals)

  filePath.write(',%d,%d' % (lowerIntervalAverages, upperIntervalAverages))

def generateHeterozygotes(numIndividuals):
  generation = []
  Aparent = Diploid(name = "A", newChr = 6)
  Bparent = Diploid(name = "B", newChr = 6)

  for i in range(numIndividuals):
    generation.append(Aparent.mate(Bparent)[0])
    
  return generation
  

