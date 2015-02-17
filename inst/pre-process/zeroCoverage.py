# Find regions that have zero coverage
# Copyright (C) 2015 Harold Pimentel
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.


import array
import math
import os
import pysam
import sys

from datetime import datetime, date, time

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def print_log(aString):
    print >> sys.stderr, "[%s] %s" % (right_now(), aString)
    return

def status_bar(cur, total, barLength):
    percent = float(cur) / float(total)
    numStars = int(math.floor(percent * barLength))
    sys.stdout.write("\r[")
    for i in xrange(numStars):
        sys.stdout.write("*")
    for i in xrange(barLength - numStars):
        sys.stdout.write(" ")
    sys.stdout.write("]")
    sys.stdout.write(" " + str(cur) + " / " + str(total))
    sys.stdout.flush()
    return

def unique(seq):
   # Not order preserving
   keys = {}
   for e in seq:
       keys[e] = 1
   return keys.keys()

# make dictionary of all transcripts. when hashed, get array ength of trans
# for all reads, find the trans it maps to then add 1 to it the count
def readXprs(fileName, colName):
    """
    fileName - a filename that points to a "results.xprs" file
    colName - a column name (i.e. fpkm, eff_counts, etc...)
    Returns: a dictionary that contains the fpkm
    """
    fileHandle = open(fileName, "r")
    firstLine = fileHandle.readline()
    firstLine = firstLine.split()

    # Find the column index
    whichCol = -1
    for colIdx in xrange(len(firstLine)):
        if (firstLine[colIdx] == colName):
            whichCol = colIdx
            break

    if whichCol == -1:
        print_log("Error: " + colName + " is not a valid column name")
        return

    xprsDict = {}
    lower_thresh = 5.0/1000.0
    for line in fileHandle:
        line = line.split()
        if float(line[4]) != 0:
            #xprsDict[line[1]] = float(line[whichCol])

            # make sure the rate is bigger than the prior
            rate = float(line[6]) / float(line[2])
            if rate > lower_thresh:
                # compute the rate: expected counts per effective base
                # use the effective length for the rate used in the computations
                # print >> sys.stderr, line[1], rate
                xprsDict[line[1]] = float(line[6]) / float(line[3])

    fileHandle.close()

    return xprsDict

# @profile
def computeZeroRegions(trans, rate):
    # utrans = unique(trans)
    utrans = trans.keys()
    utrans.sort()
    # utrans.sort()
    ps = []
    ranges = []
    L = utrans[0] - 1
    # if L >= 1:
    # print trans
    if L > 1:
        # ps.append( math.exp(-rate*L) )
        ps.append( -rate*L ) # taking the log
        ranges.append( (1, utrans[0] - 1) )
    for i in xrange(len(utrans) - 1):
        L = utrans[i + 1] - utrans[i] - 1
        # if L >= 1:
        if L > 1:
            # ps.append( math.exp(-rate*L) )
            ps.append( -rate*L ) # taking the log pval
            ranges.append( (utrans[i] + 1, utrans[i + 1] - 1) )

    if len(ps) == 0:
        return None
    return (ranges, ps)

def computeStartingPos(xprsResults, bamFile):
    # create an empty list of lists
    trans = [[] for x in xrange(len(bamFile.lengths))]
    try:
        while True:
            read = bamFile.next()
            if read.tid != -1:
                # trans[read.tid].append(read.qstart)
                trans[read.tid].append(read.pos)
    except StopIteration:
        pass

    return trans

# faster than its counterpart
def computeStartingPos2(xprsResults, bamFile):
    # create an empty list of lists
    trans = [{} for x in xrange(len(bamFile.lengths))]
    try:
        while True:
            read = bamFile.next()
            if read.tid != -1:
                # read.pos += 1 # make it 1 based
                startPos = read.pos + 1 # make it 1 based
                if read.is_reverse:
                    startPos = startPos + read.rlen - 1
                if startPos in trans[read.tid]:
                    trans[read.tid][startPos] += 1
                else:
                    trans[read.tid][startPos] = 1
    except StopIteration:
        pass

    return trans

# @profile
def computeTests(trans, xprsResults, bamFile):
    ntrans = len(trans)
    tests = [None] * ntrans
    for i in xrange(ntrans):
        t = trans[i]
        curTest = None
        if i % 50 == 0:
            status_bar(i, ntrans, 50)
        if bamFile.getrname(i) in xprsResults:
            curTest = computeZeroRegions(t, xprsResults[ bamFile.getrname(i) ])
        tests[i] = curTest
    print # clear the screen
    return tests

# @profile
def printRegions(tests, outHandle, bamFile):
    outHandle.write("reference\tstart\tend\tpvalue\n")
    for tid in xrange(len(tests)):
        # for each test per transcript
        if tests[tid] is None or len(tests[tid]) != 2:
            continue
        status_bar(tid, len(tests), 50)
        target = bamFile.getrname(tid)
        for (range, p) in zip(tests[tid][0], tests[tid][1]):
            outLine = target + "\t" + str(range[0]) + "\t" + str(range[1]) + "\t" + str(p) + "\n"
            outHandle.write(outLine)
    print

def main():
    """docstring for main"""
    # read xprs data -- returns dict with fpkm
    xprsFName = sys.argv[1]
    print_log("Reading file " + xprsFName)
    xprsResults = readXprs(xprsFName, "fpkm")

    # read bam header -- returns hash with count array
    # XXX: Only BAM currently
    samFName = sys.argv[2]
    print_log("Opening BAM file " + samFName)
    samFile = pysam.Samfile(samFName, "rb")

    print_log("Compute the number of reads starting at each position in the transcriptome")
    trans = computeStartingPos2(xprsResults, samFile)

    print_log("Computing Pr[zero region]")
    tests = computeTests(trans, xprsResults, samFile)

    outHandle = open(sys.argv[3], "w")
    print_log("Printing...")
    printRegions( tests, outHandle, samFile )
    outHandle.close()

    samFile.close()

if __name__ == '__main__':
    main()

# python ~/zeroCoverage/zeroCoverage.py ~/er/human/ec-jl-794/Sample_Orthochromatic/xprs/results.xprs ~/er/human/ec-jl-794/Sample_Orthochromatic/align.bam ~/zeroCoverage/tests.out
# python -m cProfile ~/zeroCoverage/zeroCoverage.py ~/er/human/ec-jl-794/Sample_Orthochromatic/xprs/results.xprs ~/er/human/ec-jl-794/Sample_Orthochromatic/align.bam ~/zeroCoverage/tests.out
# kernprof.py -l -o zcHash.lprof ~/zeroCoverage/zeroCoverage.py ~/er/human/ec-jl-794/Sample_Orthochromatic/xprs/results.xprs ~/er/human/ec-jl-794/Sample_Orthochromatic/align.bam ~/zeroCoverage/tests.out.hash
