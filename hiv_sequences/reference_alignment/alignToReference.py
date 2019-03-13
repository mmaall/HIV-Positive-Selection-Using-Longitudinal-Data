#!/bin/bash/python

import sys
import argparse
import getopt
import os
import glob
import shlex
import subprocess

#Returns a list of triplets a,b,c
#a: unique id
#b: sequence header
#c: sequence

verbose=False

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
        fileH.close()
        yield header,sequence

def trimEndOfSequence(sequence, char):
    for currChar in reversed(sequence):
        if(currChar== char):
            sequence=sequence[:-1]
        else: break
    return sequence

def fastaToSequenceList(fileName):

    file=None
    try:
    	file=open(fileName, "r")
    except Exception: 
    	raise Exception("Open: "+fileName+": Unable to open")  

    fileSequenceList= []

    fastaFileName=''
    sequence= ''
    header=''
    #Iterate through the tail
    for line in file:
        if line[0]==">":
            if sequence != '' and fastaFileName != '' :
                #We need to add the file name and sequence to the 
                
                #must trim  dashs (-) from the end
                sequence= trimEndOfSequence(sequence, '-')


                #file sequence list
                fileSequenceList.append((fastaFileName,header, sequence))
                sequence=''
                fastaFileName=''
                header=line
            else:
      	        header=line

            #We have reached a new sequence let's read
            firstUnderScore= False
            for char in line:
                if char == '>' : continue #ignore first >
                if char == '_' and firstUnderScore:
                    break
                    #Unique identifier read in
                elif char == '_' and not firstUnderScore :
                    #Found first underscore
                    firstUnderScore= True
                #add the current
                fastaFileName+= char

        else:
            #add current line to sequence
            sequence+= line.rstrip('\n')

    #add in last trailing file
    fileSequenceList.append((fastaFileName,header ,sequence))
    file.close()
    return fileSequenceList



def alignToReference(globalAlignmentFile, fastaDir):
    globalSequences=[]
    fastaTail= ".fasta_linsi.fasta"

    try:
        reader= FastaReader(globalAlignmentFile)
        for header,seq in reader.readFasta():
            globalSequences.append((header,seq))
    except Exception:
        raise Exception("verifyAligned: "+globalAlignmentFile+": Unable to open")
    
    #print(globalSequences[0][0])
    #print(globalSequences[0][1])
    fastaFileNameList = glob.glob(fastaDir+'*'+fastaTail)
    
    sequenceList = []
    
    for fileName in fastaFileNameList:
        reader = FastaReader(fileName)
        allSequences = []
        for header, seq in reader.readFasta():
            allSequences.append((header,seq))
        
        finalSeq= allSequences[-1]
        file= open(fileName, "w")
        file.write(">"+globalSequences[0][0]+"\n")
        file.write(globalSequences[0][1]+"\n")
        file.write(">"+finalSeq[0]+"\n")
        file.write(finalSeq[1]+"\n")
        file.close()
        subprocess.run(['muscle','-in', fileName, "-out", fastaDir+fileName])


def trim(fastaDir, size ):
    fastaTail= ".fasta_linsi.fasta"
    fastaFileNameList = glob.glob(fastaDir+'*'+fastaTail)
    for fileName in fastaFileNameList: 
        reader = FastaReader(fileName)
        seqList = []
        for header,sequence in reader.readFasta():
            seqList.append((header,sequence))
        file = open(fileName, "w")
        for header, sequence in seqList:
            newSeq= sequence [size:]
            file.write(">"+header+"\n")
            file.write(newSeq+"\n")

        file.close()



def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument( '-d', dest='fastaDir', help='Fasta File Directory', required=True)
    p.add_argument('-a', dest='alignmentFile', help='Global alignment file', required=True)
    p.add_argument('-v', '--verbose', help= 'Verbose printing',action='store_true')
    args = p.parse_args()

    if args.verbose :
    	global verbose
    	verbose=True

    if args.fastaDir is None:
        p.print_help()

    if args.alignmentFile is None:
        p.print_help()

    alignToReference(args.alignmentFile, args.fastaDir)
    size= 2549
    trim(args.fastaDir, size)





if __name__ == '__main__':
	main(sys.argv[1:])