#!/bin/bash/python


import sys
import argparse
import getopt
import os
from Bio import AlignIO

#Used to align fasta files
#Uses global allignment file provided to re allign all the patient
#patient files. Allignment done using BioPython Libarary 


def alignFiles(alignment,fastaDir):
    #Holds the tail of the fasta file name
    fastaTail= ".fasta_linsi.fasta"
    alignmentFile= open(alignment, "r")

    #Contains a triplet (a,b,c)
    #a: the fasta file name
    #b: the header for the sequence
    #c: the sequence itself
    fileSequenceList= []

    fastaFileName=''
    sequence= ''
    header=''
    #Iterate through the tail
    for line in alignmentFile:


        if line[0]==">":
            if sequence != '' and fastaFileName != '' :
                #We need to add the file name and sequence to the 
                #file sequence list

                fileSequenceList.append((fastaFileName,header, sequence))
                sequence=''
                fastaFileName=''
                header=line

            #We have reached a new fasta file let's read
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
            #Add standard tail to fastaFileName
            fastaFileName+= fastaTail



        else:
            #add current line to sequence
            sequence+= line.rstrip('\n')

    #add in last trailing file
    fileSequenceList.append((fastaFileName,header ,sequence))

    """
    for a,b,c in fileSequenceList:
        print (a,b)
        file = open(fastaDir+a, "r")
        print ('Opened '+ fastaDir+ a)
        file.close()
        print ('Closed '+ fastaDir+a)
        print()
        print()
    """

    #create folder holding global alignment
    

    alignDir= "./global_alignment/"

    try:
        os.mkdir(alignDir)
    except FileExistsError:
        print ("Mkdir: Directory " +alignDir+ " already exists")

    for file, head, seq in fileSequenceList:
        #Iterate through all members of the file sequence list
        oldAlignment= None
        newAlignment= None
        try: 
            #open old alignment file
            oldAlignment= open(fastaDir+file, "r")
        except Exception:
            print ("Open: Unable to open "+fastaDir+file)
            print ("File not aligned")
            continue

        try:
            #create new alignment file
            newAlignment= open(alignDir+file, "w")
        except Exception:
            print("Open: Unable to create" +alignDir+file)
            print("File not aligned")

        #Copy in globaly aligned data
        #write in header
        newAlignment.write(head)
        #write in sequences, adding a newline every 60 characters
        count= 0
        for char in seq:
            newAlignment.write(char)
            if count==58 :
                newAlignment.write("\n")
                count=0
            else:
                count+=1
        newAlignment.write("\n\n")

        #write in old alignment
        numSequences=0
        for line in oldAlignment:
            if line[0] == ">" :
                numSequences+=1

            if numSequences > 1:
                newAlignment.write(line)

        newAlignment.close()
        oldAlignment.close()


        #align =AlignIO.read(open(alignDir+file, "r"), "fasta")
        #print (align)












def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument( '-d', dest='fastaDir', help='Fasta File Directory', required=True)
    p.add_argument('-a', dest='alignmentFile', help='Global alignment file', required=True)
    args = p.parse_args()


    if args.fastaDir is None:
        p.print_help()

    if args.alignmentFile is None:
        p.print_help()


    alignFiles (args.alignmentFile, args.fastaDir)





if __name__ == '__main__':
	main(sys.argv[1:])