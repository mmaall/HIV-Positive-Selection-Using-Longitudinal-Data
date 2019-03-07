#!/bin/bash/python

import sys
import argparse
import getopt
import os

#Returns a list of triplets a,b,c
#a: unique id
#b: sequence header
#c: sequence

verbose=False


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


#Verifies whether the patient files provided are aligned with the
#global alignment file
#
#Args
#globalAlignmentFile: file with globaly aligned sequences
#fastaDir: Directory holding all patients fasta files
#			Must be in ./<path>/ notation. Must end in /
def verifyAligned(globalAlignmentFile, fastaDir):

    fastaTail= ".fasta_linsi.fasta"
    globalSequences= None
    try:
    	globalSequences = fastaToSequenceList(globalAlignmentFile)
    except Exception:
    	raise Exception("verifyAligned: "+globalAlignmentFile+": Unable to open")

    fileNotFoundCount= 0
    failureCount=0 
    for globalID, globalHeader, globalSeq in globalSequences:

        localPath= fastaDir+globalID+fastaTail
        localSequences=None
        try:
        	localSequences= fastaToSequenceList(localPath)
        except Exception:
        	print("File: "+localPath+": File Not Found")
        	print("Unable to analyze "+globalID)
        	fileNotFoundCount+=1


        for localID, localHeader, localSeq in localSequences:
            if globalHeader == localHeader:
                if globalSeq == localSeq:
                	if verbose:
                   	    print("Sequence " + globalID +" matches!")

                else:
                    failureCount+=1
                    print("FAILURE: Sequence: " + globalID + " does not match!")
                    if verbose:
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

    print("Number of unmatched sequences: {}".format(failureCount))	
    print("Number of sequences not evaluated: {}".format(fileNotFoundCount))


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

    verifyAligned (args.alignmentFile, args.fastaDir)





if __name__ == '__main__':
	main(sys.argv[1:])