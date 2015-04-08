#!/usr/bin/env python
# encoding: utf-8
"""
Chromosomes.py

Created by Joshua Shapiro on 2008-08-11.
"""
from __future__ import division
from numpy import *
from WormUtils import *
import itertools
import operator

def generateBreaksPoisson(cM = 200):
    breaks = random.uniform(size = random.poisson(cM/100.0))
    breaks.sort()
    return breaks

class Chromosome(object):
    """Chromosome object which contains information on parentage of segments"""
    def __init__(self,  cM=200, name=None, segments=None, newParent=None, interference = "absent"):
        super(Chromosome, self).__init__()
        self.name = name
        self.cM = cM
        self.segments = segments
        #segments are lists of tuples where the first position is the start location, the second is the parent of origin
        self.interference = interference
        
        if newParent != None:
            if self.segments != None:
                raise ValueError, "Can't set segment information in a new parent"
            self.segments = [(0,newParent)]
        
        if not (self.interference in ["complete", "absent"]):
          raise ValueError, "Interference must be one of 'complete' or 'absent'."
        
    def __eq__(self, other):
      if  isinstance(other, Chromosome):
        if self.name != other.name:
          return False
        if len(self.segments) != len(other.segments):
          return False
        for s1, s2 in itertools.izip(self.segments, other.segments):
          if s1 != s2: return False
        return True  
      else: return NotImplemented
    
    def __ne__(self, other):
      equal_result = self.__eq__(other)
      if (equal_result is not NotImplemented):
          return not equal_result
      return NotImplemented
    
    @staticmethod 
    def getGeneticDistance(physLoc, number):
        """Gets the centimorgan location from a physical location in kilobases"""
        
        check_physLoc(physLoc, number)
        
        chromosome_breaks = [(527000, 0), (3331000, 3.43), (7182000, 1.34), (3835000, 6.78), (197000, 0),
                             (306000, 0), (4573000, 4.92), (7141000, 1.33), (2589000, 8.47), (670000, 0),
                             (494000, 0), (3228000, 7.83), (6618000, 1.17), (2877000, 7.24), (567000, 0),
                             (720000, 0), (3176000, 7.65), (9074000, 1.05), (3742000, 3.64), (782000, 0),
                             (643000, 0), (5254000, 3.22), (10653000, 1.32), (3787000, 5.47), (583000, 0),
                             (572000, 0), (5565000, 3.81), (6343000, 1.70), (3937000, 5.14), (1302000, 0)]
        genLoc = 0
        i = 0
        kb = physLoc
        
        while kb > 0:
            kb -= chromosome_breaks[(number * 5) + i][0]
            seg = chromosome_breaks[(number * 5) + i][0]
           
            if kb < 0:
                seg += kb
                
            genLoc += (seg/1000000) * chromosome_breaks[(number * 5) + i][1]
            i+=1
        return genLoc
    
    @staticmethod
    def getLoc(physLoc, number):
        genLoc = Chromosome.getGeneticDistance(physLoc, number)
        
        return genLoc / cM_max[number]
    
    @staticmethod
    def getPhysDistance(cM, number):
        """Gets the physical location in kilobases from the given centimorgan location"""
        kb_shifts_left = [527000, 306000, 494000, 720000, 643000, 572000]
        kb_shifts_right = [197000, 670000, 567000, 782000, 583000, 1302000]
        
        if cM > cM_max[number] or cM < 0:
            raise ValueError, "The centimorgan distance must be within the range of the chromosome"
        
        cM_breaks = [(11.42533, .29154), (9.62388, .74626), (26.0013, .14749),
                     (22.49916, .20325), (9.49753, .75187), (21.92883, .11806),
                     (25.27524, .12771), (7.74306, .8547), (20.82948, .13812),
                     (24.2964, .13072), (9.5277, .95238), (13.62088, .27472),
                     (16.91788, .31055), (14.06196, .75757), (20.71489, .18281),
                     (21.20265, .26246), (10.7831, .58823), (20.23618, .19455)]
        physLoc = kb_shifts_left[number]
        i = 0
        gen = cM
        while gen > 0.001:
            gen -= cM_breaks[(number * 3) + i][0]
            seg = cM_breaks[(number * 3) + i][0]
        
            if gen < 0:
                seg += gen
        
            physLoc += (seg * cM_breaks[(number * 3) + i][1])*1000000
            i += 1
        
        if cM == cM_max[number]:
            physLoc += kb_shifts_right[number]
            
        return physLoc
    
    @staticmethod
    def getPhysDistanceFromLoc(loc, number):
        """Gets a physical location from a [0,1] location"""
        if loc < 0.0 or loc > 1.0:
            raise ValueError, "The location must be within the range [0,1]"
        
        return Chromosome.getPhysDistance(loc * cM_max[number], number)
    
    def getParentAtLocation(self, loc):
        """gets the Parental Identity for a chromosomal location"""
        if loc < 0.0 or loc > 1.01:
            raise ValueError, "Location must be in range [0,1]"
        if loc == 0.0:
            return self.segments[0][1]
            
        i = 0
        while i < len(self.segments) and self.segments[i][0] < loc:
            i+=1

        return self.segments[i - 1][1]
    
    def getParentAtLocations(self, locs):
        """gets the Parental Identity for a list of chromosomal locations"""
        parents = [''] * len(locs)
        for i in range(len(locs)):
            parents[i] = getParentAtLocation(locs[i])
        
        return parents
    
    def getParentAtMapLoc(self, mapLoc):
      """gets the Parental identity for a given cM position"""
      if mapLoc < 0 or mapLoc > self.cM:
        raise ValueError, "Map location must be withing the range of the chromosome."
      loc = cM/float(self.cM)
      return self.getParentAtLocation(loc)
    
    #TODO(zifanxiang): Refactor using getParentAtMapLoc
    def getParentAtMapLocs(self, mapLocs):
        """gets the Parental Identity for a list of chromosomal locations"""
        parents = [''] * len(mapLocs)
        locs = [ml/float(self.cM) for ml in mapLocs]
        order = [i for _,i in sorted(itertools.izip( locs, range(len(locs)) ))]
        if locs[order[0]] < 0 or locs[order[-1]] > 1 :
            raise ValueError, "Locations must be in range [0,1]"
        i = 0
        for n in order:
          while i < len(self.segments) and self.segments[i][0] < locs[n] :
            i+=1
          parents[n] = self.segments[i-1][1]
        return parents
    
    def recombine(self, mate, interference = None):
        if self.name != mate.name:
            raise ValueError, "Chromosome names are not the same; can't recombine between them." 
        if self == mate:
          #shortcut: any recombinants would be identical anyway
          return (Chromosome(name = self.name, cM = self.cM,  segments = self.segments), 
                  Chromosome(name = self.name, cM = self.cM,  segments = self.segments))
        
        if interference == None:
          interference = self.interference
        
        segments1 = list(self.segments)
        segments2 = list(mate.segments)
        brokenSegments1 = list()
        brokenSegments2 = list()
        
        if interference == "absent":
          crossOvers = generateBreaksPoisson(self.cM)
        elif interference == "complete":
          crossOvers = random.uniform(size=1)
        else:
          raise ValueError, "Interference setting must be one of 'absent' or 'complete'."
        
        for crossOver in crossOvers:
            tempSeg = list()
            #move segments to temp list until you get past the crossover
            while len(segments1) > 0 and segments1[0][0] < crossOver:
                tempSeg.append(segments1.pop(0))
            #add in new start segment
            segments1.insert(0, (crossOver, tempSeg[-1][1]))
            brokenSegments1.append(tempSeg)
            
            tempSeg=[]
            while len(segments2) > 0 and segments2[0][0] < crossOver:
                tempSeg.append(segments2.pop(0))
            segments2.insert(0, (crossOver, tempSeg[-1][1]))
            brokenSegments2.append(tempSeg)
        #all breakpoints done, now just add the remainder of the segments to our lists
        brokenSegments1.append(segments1)
        brokenSegments2.append(segments2)
        
        #combine breakpoint lists, alternating parental lists
        chr1 = list()
        chr2 = list()
        for index, bs1 in enumerate(brokenSegments1):
            if index % 2 == 0 : 
                chr1 += bs1 
                chr2 += brokenSegments2[index]
            else:
                chr1 += brokenSegments2[index]
                chr2 += bs1
        
        #remove redundant segments
        chr1 = [x for i,x in enumerate(chr1) if (i == 0 or chr1[i][1] != chr1[i-1][1])]
        chr2 = [x for i,x in enumerate(chr2) if (i == 0 or chr2[i][1] != chr2[i-1][1])]
        
        
        if random.binomial(1,0.5): #randomly order xover products
            return (Chromosome(name = self.name, cM = self.cM,  segments = chr1), 
                    Chromosome(name = self.name, cM = self.cM,  segments = chr2))
        else:
            return (Chromosome(name = self.name, cM = self.cM,  segments = chr2), 
                    Chromosome(name = self.name, cM = self.cM,  segments = chr1))            

    #Returns the percentage of the chromosome that has the allele of a particular parent
    def getPercentageOfParent(self, parentName):
        seg = list(self.segments)
        i = 1
        percent = 0
        
        while (i < len(seg)):
            curParent = seg[i - 1][1]
            if curParent == parentName:
                percent += seg[i][0] - seg[i - 1][0]
            
            i += 1

        if seg[i - 1][1] == parentName:
            percent += 1.0 - seg[i - 1][0]
        
        return percent
    
    #Returns the left and right bounds of the allele at the selected genetic location
    def physicalLocsOfInterval(self, loc, chromNumber):
        seg = list(self.segments)
        i = 0
        
        while (i < len(seg) and seg[i][0] < loc):
            i += 1
    
        lowerBound = 0
        upperBound = 1
        
        if (i < len(seg)):
            upperBound = seg[i][0]
    
        if (i > 0):
            lowerBound = seg[i - 1][0]
        elif (len(seg) != 1):
            upperBound = seg[i + 1][0]
        else:
	        upperBound = 1
        
        return [int(Chromosome.getPhysDistanceFromLoc(lowerBound, chromNumber)), int(Chromosome.getPhysDistanceFromLoc(upperBound, chromNumber))]

    def getLengthGeneticSizes(self, chromNumber):
        seg = list(self.segments)
        lengthGeneticSizes = []

        for i in range(1, len(seg)):
            locSize = seg[i][0] - seg[i - 1][0]
            geneticSize = locSize * cM_max[chromNumber]
            lengthGeneticSizes.append(locSize)

        locSize = 1.0 - seg[len(seg) - 1][0]
        lengthGeneticSizes.append(locSize * cM_max[chromNumber])

        return lengthGeneticSizes
        
if __name__ == '__main__':
    #x = Chromosome.getLoc(15072, 0)
    #print(x)
    #y = Chromosome.getPhysDistanceFromLoc(x, 0)
    #print(y)
    #z = Chromosome.getPhysDistanceFromLoc(.5, 1)
    #print(z)
    u = Chromosome.getPhysDistance(47.0507, 0)
    print(u)
    #a = Chromosome(newParent = "Blue")
    #b = Chromosome(segments = [(0,1)])
    #newChrs = a.recombine(b)
    #print newChrs[0].segments
    #newChrs = newChrs[0].recombine(newChrs[1])
    #k = newChrs[0].physicalLocsOfInterval(0, 0)
    #k2 = newChrs[0].physicalLocsOfInterval(.6, 0)
    #k3 = newChrs[0].physicalLocsOfInterval(.8, 0)
    #print "Chromosome Set 1: ", (newChrs[0].segments)
    #print "Chromosome Set 2: ", (newChrs[1].segments)
    #print "lower: %d upper: %d\n" % (k[0], k[1])
    #print "lower: %d upper: %d\n" % (k2[0], k2[1])
    #print "lower: %d upper: %d\n" % (k3[0], k3[1])
    #print newChrs[0].segments
    #print newChrs[0].getPercentageOfParent("Blue")
    #print (list(newChrs[0].segments))
    #print newChrs[0].segments
    #print newChrs[0].getParentAtLocation(0.5)