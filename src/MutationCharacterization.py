#!/usr/bin/env python3
"""
Characterization of mutations between sequences
"""

from enum import Enum
#from FastaReader import FastaReader

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

        # Key = Position, value = [codon, Yt, Yv, St, Sv]
        self.mutCharDict = {}

    def findMutations(self):
        mutType = ""
        codon = 0 # to start off
        for pos in range(0, len(self.seqt0), 3):
            t0codon = self.seqt0[pos:pos+3]
            tfcodon = self.seqtf[pos:pos+3]
            codon += 1 # we want to know which codon it occurs in
            if t0codon != tfcodon:
                if t0codon in self.dnaCodonTable.keys() and tfcodon in self.dnaCodonTable.keys():
                    if self.dnaCodonTable[t0codon] != self.dnaCodonTable[tfcodon]: #nonsyn

                        for i in range(3):
                            if t0codon[i] != tfcodon[i]: 
                                if self.baseStructure[t0codon[i]] != self.baseStructure[tfcodon[i]]:
                                    # nonsyn and transversion
                                    self.mutCharDict.update({(pos+i+1):[codon, 0, 1, 0, 0]})
                                else:
                                    # nonsyn and transition
                                    self.mutCharDict.update({(pos+i+1):[codon, 1, 0, 0, 0]})

                    else: # synonymous
                        for i in range(3):
                            if t0codon[i] != tfcodon[i]: 
                                if self.baseStructure[t0codon[i]] != self.baseStructure[tfcodon[i]]:
                                    # syn and transversion
                                    self.mutCharDict.update({(pos+i+1):[codon, 0, 0, 0, 1]})
                                else:
                                    # syn and transition
                                    self.mutCharDict.update({(pos+i+1):[codon, 0, 0, 1, 0]})

        return self.mutCharDict

class PatientProfile :
    def __init__ (self, header, mutDatabase):
        thisHeader = header.split('_')
        self.drugs = []
        # records all drugs taken by patient
        for i in range(6, len(thisHeader)-1):
            if len(thisHeader[i]) > 1: 
                self.drugs.append(thisHeader[i])

        
        self.patientID = thisHeader
        self.mutations = mutDatabase
        self.transverseCount = 0
        self.transitionCount = 0
        self.synCount = 0
        self.nonsynCount = 0
        
        # Saves count for t, v to be used in calculating t,v frequencies
        for mutInfo in self.mutations.values():
            if mutInfo[2] == 1 or mutInfo[4]:
                self.transverseCount += 1
            elif mutInfo[1] == 1 or mutInfo[3] == 1:
                self.transitionCount += 1

        for mutInfo in self.mutations.values():
            if mutInfo[1] == 1 or mutInfo[2] == 1:
                self.nonsynCount += 1
            elif mutInfo[3] == 0 or mutInfo[4] == 1:
                self.synCount += 1


import sys
class FastaReader :
    
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


def main():

    myReader = FastaReader()
    seqDict = []
    patientInfo = {}
    for header, sequence in myReader.readFasta():
        seqDict.append([header, sequence])

    myMutChar = mutCharStorage(seqDict[0][1], seqDict[1][1])
    totalMutations = myMutChar.findMutations()
    thisPatient = PatientProfile(seqDict[1][0], totalMutations)
    print(thisPatient.patientID)
    print("drugs taken: "+str(thisPatient.drugs))
    print("v: "+str(thisPatient.transverseCount))
    print("t: "+str(thisPatient.transitionCount))
    print("nonsyn: "+str(thisPatient.nonsynCount))
    print("syn: "+str(thisPatient.synCount))
    for key in totalMutations:
        print(str(key)+": "+str(totalMutations[key]))

if __name__ == '__main__':
    main()
