#!/bin/bash/python

import sys
from enum import Enum
from FastaReader import FastaReader

#Creates an enum for mutation types
#Contains transitions and transversions
class transTranv(Enum):
    transition= 1
    transversion= 2

#Creates an enum for types of bases
#Contains both purines and pyrimidines
class mutationType(Enum):
    syn= 1
    nonsyn= 2

class baseType(Enum):
    pyr= 1
    pur=2


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
    'H': ['A','C','T'],
	'N': ['A','C','G','T']

}

baseStructure = {'C': baseType.pyr, 'T': baseType.pyr, 'G': baseType.pur, 'A': baseType.pur}
   
def findMutations(seqt0, seqtf):
    mutType = ""
    # Key = Position, value = [t0 base, tf base, transition/transversion, mutType]
    mutCharDict= {}
    for pos in range(0, len(seqt0), 3):
        t0codon = seqt0[pos:pos+3]
        tfcodon = seqtf[pos:pos+3]
        foundAmbiguousBase= False
        for pos1, pos2 in zip(t0codon, tfcodon):
            if pos1 in nucleicAcidCodeTable or pos2 in nucleicAcidCodeTable:
                foundAmbiguousBase=True
                break

        if foundAmbiguousBase:
            continue

        if t0codon != tfcodon:
            if "-" in t0codon or "-" in tfcodon:
                continue
            if dnaCodonTable[t0codon] == dnaCodonTable[tfcodon]:
                mutType = mutationType.syn
            else:
                mutType = mutationType.nonsyn
            for i in range(3):
                if t0codon[i] != tfcodon[i]:
                    if baseStructure[t0codon[i]] != baseStructure[tfcodon[i]]:
                        mutCharDict.update({(pos+i+1):[t0codon[i], tfcodon[i], transTranv.transversion, mutType]})
                    else:
                        mutCharDict.update({(pos+i+1):[t0codon[i], tfcodon[i], transTranv.transition, mutType]})

    return mutCharDict


def findAllPossibleMutations(seqInit):
    base= ['A', 'T', 'G','C']
    possibleMutationList= [None] * (len(seqInit)+1) 
    mutType = ""
    # Key = Position, value = [t0 base, tf base, transition/transversion, mutType]
    for pos in range(0, len(seqInit), 3):
        t0codon = seqInit[pos:pos+3]
        foundAmbiguousBase= False
        for pos1 in t0codon:
            if pos1 in nucleicAcidCodeTable or pos1 == '-':
                foundAmbiguousBase=True
                break


        if foundAmbiguousBase:
            continue


        for changedPos, currBase in enumerate(t0codon):
            #print("findAllPossibleMutations: position: "+str(changedPos+pos))
            synTransition= 0
            nonsynTransition= 0
            synTransversion= 0
            nonsynTransversion= 0    
            codonList = list(t0codon)
            for testBase in base:
                if testBase == currBase:
                    #Base to be tested is same as original base
                    continue
                isSyn=False
                isTrans=False
                codonList[changedPos]=testBase
                testCodon= ''.join(codonList)
                if dnaCodonTable[t0codon] == dnaCodonTable[testCodon]:
                    isSyn=True
                if baseStructure[t0codon[changedPos]]==baseStructure[testCodon[changedPos]]:
                    isTrans=True

                if isSyn and isTrans:
                    synTransition+=1
                elif isSyn and not isTrans:
                    synTransversion+=1
                elif not isSyn and isTrans:
                    nonsynTransition+=1
                elif not isSyn and not isTrans:
                    nonsynTransversion+=1

            possibleMutationList[pos+changedPos+1]= (synTransition,nonsynTransition, synTransversion, nonsynTransversion)
    return possibleMutationList


#Finds whether there is a mutation or not at a specific base position
#in a codon. 
#BasePosition: The position of the base in a codon(0,1,2)
#codonInit: Three letter string of the initial codon
#codonFinal: Three letter string of the final codon





"""
Patient Class
Holds patient data from a single fasta patient file.

Instance variables
fname: File Name
uniqueID: Unique Identifier provided by the file name
seqt0: A string that contains the initial sequence from the patient
seqtf: A string that contains the final sequence from the patient
mutCharDict: A dictionary that holds all the mutations from t0 to tf
    Key: Position the position of the mutation (1 to len(seqt0))
    Value: [t0 base, tf base, transition/transversion, mutType]
drugsGiven: A list of the drugs administered to this patient
possibleMutation: A list of len(seqt0) that contains a tuple of the counts of all possible mutations
    tuple (synonymous transitions, nonsynonymous transitions, synonymous transversions, nonsynonymous transversions)
"""
class Patient :
    
    def __init__ (self, fname=None,uniqueID= None, seqt0=None,
                     seqtf= None, mutCharDict= None, drugsGiven= None, possibleMutations=None):
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

        if possibleMutations== None:
            possibleMutations= []
        else:
            self.possibleMutations= possibleMutations

        self.transitionCount= 0
        self.transversionCount= 0
        self.effectiveLength= 0

    def findMutations(self):
        self.transitionsCount= 0
        self.transversionCount= 0
        self.effectiveLength=0
        mutType = ""
        # Key = Position, value = [t0 base, tf base, transition/transversion, mutType]
        self.mutCharDict= {}
        for pos in range(0, len(self.seqt0), 3):
            t0codon = self.seqt0[pos:pos+3]
            tfcodon = self.seqtf[pos:pos+3]
            foundAmbiguousBase= False
            currLength= 0
            for pos1, pos2 in zip(t0codon, tfcodon):
                if pos1 in nucleicAcidCodeTable or pos2 in nucleicAcidCodeTable:
                    foundAmbiguousBase=True
                    break
                currLength+=1
            self.effectiveLength+=currLength

            if foundAmbiguousBase:
                continue

            if t0codon != tfcodon:
                if "-" in t0codon or "-" in tfcodon:
                    continue
                if dnaCodonTable[t0codon] == dnaCodonTable[tfcodon]:
                    mutType = mutationType.syn
                else:
                    mutType = mutationType.nonsyn
                for i in range(3):
                    if t0codon[i] != tfcodon[i]:
                        if baseStructure[t0codon[i]] != baseStructure[tfcodon[i]]:
                            self.mutCharDict.update({(pos+i+1):[t0codon[i], tfcodon[i], transTranv.transversion, mutType]})
                            self.transversionCount+=1
                        else:
                            self.mutCharDict.update({(pos+i+1):[t0codon[i], tfcodon[i], transTranv.transition, mutType]})
                            self.transitionCount+=1

        

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

        self.findMutations()
        self.possibleMutations= findAllPossibleMutations(self.seqt0)
        #Parse the header and put in relevant information
        finalHeader= mutationList[-1][0]
        #print(finalHeader)
        readHeader= True
        firstUnderScore= True
        builtStr=''
        readDrugs=False
        for char in header:
            #print("Char:" +char)
            #print("builtStr: "+builtStr)
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
        MAX_LINE_LENGTH= 60
        output+="\nInitial Sequence:\n"
        for char in self.seqt0:
            if count == MAX_LINE_LENGTH:
                output+="\n"
                count= 0
            output+=char
            count+=1

        output+="\nFinal Sequence:\n"
        count= 0
        for char in self.seqtf:
            if count == MAX_LINE_LENGTH:
                output+="\n"
                count= 0
            output+=char  
            count+=1

        output+="\nDrugs Administered:\n"
        output+= "\t" +str(self.drugsGiven)
        output+="\nMutations Identified:\n"
        for key in self.mutCharDict:
            output+= "\tKey: " +str(key) +"\n"
            output+= "\tValue: "+ str(self.mutCharDict[key])+"\n"

        output+="Mutation Count: \n"
        output+="\tTransversions:\t" +str(self.transversionCount)+"\n"
        output+="\tTransitions:\t" +str(self.transitionCount) +"\n"

        output+="Effective Length:\n"
        output+= "\t"+str(self.effectiveLength) +"\n"

        #This isn't really helpful stuff, so we'll just not worry about it for now
        '''
        output+="\nPossible Mutation Counts:\n"
        for cnt, possibleChanges in enumerate(self.possibleMutations):
            if(possibleChanges== None):
                output+= "\t" +str(cnt) +": No counts\n"
            else:
                output+="\t"+str(cnt)+": " + str(possibleChanges) + "\n"
        '''


        return output



def main(argv):
    patient= Patient()
    patient.inputFile("../hiv_sequences/aligned/testInput")
    print(patient)




if __name__ == '__main__':
    main(sys.argv[1:])



