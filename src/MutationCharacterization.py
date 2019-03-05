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

class mutCharStorage:

    dnaCodonTable = {
        # DNA codon table
        # T
        'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',  # TxT
        'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',  # TxC
        'TTA': 'L', 'TCA': 'S', 'TAA': '-', 'TGA': '-',  # TxA
        'TTG': 'L', 'TCG': 'S', 'TAG': '-', 'TGG': 'W',  # TxG
        # C
        'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',  # CxT
        'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',  # AxT
        'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',  # GxT
        'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
        }

    baseStructure = {'C': 'pyr', 'T': 'pyr', 'G': 'pur', 'A': 'pur'}
    ambigCodes = {'Y': ['C', 'T'],
                  'R': ['A', 'G'],
                  'W': ['A', 'T'],
                  'S': ['G', 'C'],
                  'K': ['T', 'G'],
                  'M': ['C', 'A'],
                  'D': ['A', 'G', 'T'],
                  'V': ['A', 'C', 'G'],
                  'H': ['A', 'C', 'T'], 
                  'B': ['C', 'G', 'T']
    }

    def __init__(self, t0Seq, tfSeq):
        
        self.seqt0 = t0Seq
        self.seqtf = tfSeq

        # Key = Position, value = [t0 base, tf base, transition/transversion, mutType]
        self.mutCharDict = {}

    def findMutations(self):
        mutType = ""
        for pos in range(0, len(self.seqt0), 3):
            t0codon = self.seqt0[pos:pos+3]
            tfcodon = self.seqtf[pos:pos+3]
            if t0codon != tfcodon:
                # checks if any ambiguous codes are present in either codon
                for base in ambigCodes.keys():
                    # if ambiguous code in t0codon
                    if base in t0codon:
                        for i in range(len(ambigCodes[base])):
                            testcodon = t0codon.replace(base, ambigCodes[base][i])
                            if self.dnaCodonTable[testcodon] != self.dnaCodonTable[tfcodon]:
                                mutType = "nonsyn"
                                possibleMutation = pos+testcodon.find(base)+1
                            else:
                                syn += 1 # if it is a synonymous mutation, we don't want to count that mutation position.
                                mutType = "syn"
                            if syn > 0:
                                break # breaks out of this for loop
                            else:
                                if self.baseStructure[base] != self.baseStructure[tfcodon[possibleMut]]:
                                    self.mutCharDict.update({possibleMutation:[base, tfcodon[i], "transversion", mutType})
                                else:
                                    self.mutCharDict.update({possibleMutation:[base, tfcodon[i], "transition", mutType})
                    # if amibugous code is present in tfcodon
                    elif base in tfcodon:
                        for i in range(len(ambigCodes[base]:
                                           testcodon = tfcodon.replace(base, ambigCodes[base][i])
                            if self.dnaCodonTable[testcodon] != self.dnaCodonTable[t0codon]:
                                nonsyn += 1
                                mutType = "nonsyn"
                                possibleMutation = pos+testcodon.find(base)+1
                            else:
                                syn += 1
                                mutType = "syn"
                            if syn > 0:
                                break
                            else:
                                if self.baseStructure[base] != self.baseStructure[t0codon[possibleMut]]:
                                    self.mutCharDict.update({possibleMutation:[t0codon[i], base, "transversion", mutType})
                                else:
                                    self.mutCharDict.update({possibleMutation:[t0codon[i], base, "transition", m
                # assumes there are no ambiguous bases in either codons
                if self.dnaCodonTable[t0codon] != self.dnaCodonTable[tfcodon]:
                    mutType = "nonsyn"
                else:
                    mutType = "syn"
                for i in range(3):
                    if t0codon[i] != tfcodon[i]:
                        if self.baseStructure[t0codon[i]] != self.baseStructure[tfcodon[i]]:
                            self.mutCharDict.update({(pos+i+1):[t0codon[i], tfcodon[i], "transversion", mutType]})
                        else:
                            self.mutCharDict.update({(pos+i+1):[t0codon[i], tfcodon[i], "transition", mutType]})
        return self.mutCharDict


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
