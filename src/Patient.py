#!/bin/bash/python

import sys
from enum import Enum
from FastaReader import FastaReader

#Creates an enum for mutation types
#Contains transitions and transversions
class transitionTransversion:
    transition= 1
    transversion= 2

#Creates an enum for types of bases
#Contains both purines and pyrimidines
class mutationType:
    syn= 1
    nonsyn= 2



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

#Nucleic acid codon table
#Key: Nucleic acid code
#Value: List of bases they code to 
#Does not return anything if an ATGC
nucleicAcidCodeTable = {
	'R': ['A','G'],
	'Y': ['C','T'],
	'K': ['G','T'],
	'M': ['A','C'],
	'S': ['C','G'],
	'W': ['A','T'],
	'B': ['C','G','T'],
	'D': ['A','G','T'],
	'V': ['A','C','G'],
	'N': ['A','C','G','T']
}

baseStructure = {'C': 'pyr', 'T': 'pyr', 'G': 'pur', 'A': 'pur'}
   
def findMutations(seqt0, seqtf):
    mutType = ""
    # Key = Position, value = [t0 base, tf base, transition/transversion, mutType]
    mutCharDict= {}
    for pos in range(0, len(seqt0), 3):
        t0codon = seqt0[pos:pos+3]
        tfcodon = seqtf[pos:pos+3]
        if t0codon != tfcodon:
            if "-" in t0codon or "-" in tfcodon:
                continue
            if dnaCodonTable[t0codon] != dnaCodonTable[tfcodon]:
                mutType = "nonsyn"
            else:
                mutType = "syn"
            for i in range(3):
                if t0codon[i] != tfcodon[i]:
                    if baseStructure[t0codon[i]] != baseStructure[tfcodon[i]]:
                        mutCharDict.update({(pos+i+1):[t0codon[i], tfcodon[i], "transversion", mutType]})
                    else:
                        mutCharDict.update({(pos+i+1):[t0codon[i], tfcodon[i], "transition", mutType]})
    return mutCharDict

#Finds whether there is a mutation or not at a specific base position
#in a codon. 
#BasePosition: The position of the base in a codon(0,1,2)
#codonInit: Three letter string of the initial codon
#codonFinal: Three letter string of the final codon
def mutationAtBase(basePosition, codonInit, codonFinal):
    isMuation= False


def findMutationsComplex(seqt0, seqtf):
    mutType = ""
    # Key = Position, value = [t0 base, tf base, transition/transversion, mutType]
    mutCharDict= {}
    for pos in range(0, len(seqt0), 3):
        t0codon = seqt0[pos:pos+3]
        tfcodon = seqtf[pos:pos+3]
        if t0codon != tfcodon:
            if "-" in t0codon or "-" in tfcodon:
                continue
            if dnaCodonTable[t0codon] != dnaCodonTable[tfcodon]:
                mutType = "nonsyn"
            else:
                mutType = "syn"
            for i in range(3):
                if t0codon[i] != tfcodon[i]:
                    if baseStructure[t0codon[i]] != baseStructure[tfcodon[i]]:
                        mutCharDict.update({(pos+i+1):[t0codon[i], tfcodon[i], "transversion", mutType]})
                    else:
                        mutCharDict.update({(pos+i+1):[t0codon[i], tfcodon[i], "transition", mutType]})
    return mutCharDict





class Patient :
    
    def __init__ (self, fname=None,uniqueID= None, seqt0=None , seqtf= None, mutCharDict= None, drugsGiven= None):
        '''contructor: saves attribute fname '''
        if fname == None:
            self.fname= ''
        else:
            self.fname=fname

        if uniqueID == None:
            uniqueID = ''
        else:
            self.uniqueID= uniqueID

        if seqt0 == None:
            self.seqt0= ''
        else:
            self.seqt0=seqt0

        if seqtf == None:
            self.seqtf=''
        else:
            self.seqtf=seqtf

        if mutCharDict== None:
            self.mutCharDict= {}
        else:
            self.mutCharDict= mutCharDict

        if drugsGiven == None:
            drugsGiven= []
        else:
            self.drugsGiven= drugsGiven

        
    def inputFile(self, fname):
        self.fname=fname
        self.uniqueID= ''
        self.drugsGiven=[]
        reader = FastaReader(fname)

        mutationList= []
        for header, seq in reader.readFasta():
            mutationList.append((header,seq))

        self.seqt0= mutationList[0][1]
        self.seqtf= mutationList[-1][1]
        #Shaves '>' 

        self.mutCharDict=findMutations(self.seqt0,self.seqtf)
        #Parse the header and put in relevant information
        finalHeader= mutationList[-1][0]
        print(finalHeader)
        readHeader= True
        firstUnderScore= True
        builtStr=''
        readDrugs=False
        for char in header:
            print("Char:" +char)
            print("builtStr: "+builtStr)
            if readHeader:
                if char=='_':
                    if firstUnderScore:
                        builtStr+=char
                        firstUnderScore= False
                    else:
                        readHeader= False
                        self.uniqueID= builtStr
                else:
                    builtStr+=char

            elif readDrugs:
                if char== '_':
                    self.drugsGiven.append(builtStr)
                    builtStr=''
                elif builtStr == 'None':
                    break

                else:
                    builtStr+=char

            elif builtStr== '__':
                readDrugs= True
                builtStr=''
                builtStr+=char

            elif char != '_':
                builtStr= ''
            else:
                builtStr+=char


    def __str__( self ):
        output= ''
        output+=self.fname + "\n"
        output+="Unique ID: "+self.uniqueID
        count= 0
        MAX_LINE_LENGTH= 50
        output+="\nInitial Sequence:\n"
        for char in self.seqt0:
            if count == 50:
                output+="\n"
                count= 0
            output+=char

        output+="\nFinal Sequence:\n"
        for char in self.seqtf:
            if count == 50:
                output+="\n"
                count= 0
            output+=char  

        output+="\nDrugs Administered:\n"
        output+= str(self.drugsGiven)
        output+="\nMutations Identified:\n"
        for key, values in self.mutCharDict:
            output+= "\tKey: " +key +"\n"
            output+= "\tValue: "+ values+"\n"


        return output



def main(argv):
    patient= Patient()
    patient.inputFile("../hiv_sequences/aligned/968_12525.fasta_linsi.fasta")
    print(patient)




if __name__ == '__main__':
    main(sys.argv[1:])



