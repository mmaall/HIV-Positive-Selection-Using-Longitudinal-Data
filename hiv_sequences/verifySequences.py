#!/bin/bash/python

import sys
import argparse
import getopt
import os

#Returns a list of triplets a,b,c
#a: unique id
#b: sequence header
#c: sequence


def fastaToSequenceList(fileName):

    file=open(fileName, "r")

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
                for currChar in reversed(sequence):
                	if(currChar== '-'):
                		sequence=sequence[:-1]
                	else: break


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

def verifyAligned(globalAlignmentFile, fastaDir):

    fastaTail= ".fasta_linsi.fasta"
    globalSequences = fastaToSequenceList(globalAlignmentFile)

    failureCount=0 
    for globalID, globalHeader, globalSeq in globalSequences:
        localSequences= fastaToSequenceList(fastaDir+globalID+fastaTail)
        #print("Testing "+globalID)
        for localID, localHeader, localSeq in localSequences:
            #print (globalHeader)
            #print (localHeader)
            if globalHeader == localHeader:
                if globalSeq == localSeq:
                    print("Sequence " + globalID +" matches!")

                else:
                    failureCount+=1
                    print("FAILURE: Sequence: " + globalID + " does not match!")
                    print("Global: "+globalID)
                    print(globalHeader)
                    print(globalSeq)
                    print("Local: "+localID)
                    print(localHeader)
                    print(localSeq)
                    print("\n\n\n\n")
                break

        #print (localSequences)
        #print(localSequences[0])

    


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument( '-d', dest='fastaDir', help='Fasta File Directory', required=True)
    p.add_argument('-a', dest='alignmentFile', help='Global alignment file', required=True)
    args = p.parse_args()


    if args.fastaDir is None:
        p.print_help()

    if args.alignmentFile is None:
        p.print_help()

    verifyAligned (args.alignmentFile, args.fastaDir)





if __name__ == '__main__':
	main(sys.argv[1:])