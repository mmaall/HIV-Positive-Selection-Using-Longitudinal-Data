#!/usr/bin/env python3
"""
Characterization of mutations between sequences
"""

from enum import Enum
from FastaReader import FastaReader

#Creates an enum for mutation types
#Contains transitions and transversions
class mutationType:
    transition= 1
    transversion= 2

#Creates an enum for types of bases
#Contains both purines and pyrimidines
class baseType:
    pur= 1
    pyr= 2



import sys
class FastAreader :
    
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

def main(inCL=None):

    if inCL is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(inCL)

    myReader = FastAreader()
    seqDict = []

    for header,seq in myReader.readFasta():
        seqDict.append([header, seq])

    myMutChar = mutCharStorage(seqDict[0][1], [1][1])
    totalMutations = myMutChar.findMutations()
    for key in totalMutations:
        print("key = " +totalMutations[key])

if __name__ == '__main__':
    main()