#!/bin/bash/python

import sys
import argparse
import os
import Patient

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

def fastaInfo(fileNameList, verbose):
    seqDict= []
    for fileName in fileNameList:
        reader= FastAreader(fileName)
        for header,seq in reader.readFasta():
            seqDict.append((header,seq))

    baseCount= 0
    rCount= 0 #Number of unkown purines found
    yCount= 0 #Number of unkown pyrmidines
    otherCount= 0 #Number of other bases found
    numSequences= len(seqDict)
    normalBases= {'A','T','G','C','-'}
    for header,sequence in seqDict:
        for char in sequence:
            #print(type(char))
            baseCount+=1
            if char =='R':
                rCount+=1
            elif char == 'Y':
                yCount+=1
            elif char not in normalBases:
                #print("Found abnormal: "+char)
                otherCount+=1 

    print("Number of unkown purines: "+str(rCount))
    print("Number of unkown pyrimidines: "+str(yCount))
    print("Number of other non-AGTC codes: "+ str(otherCount))
    print("Total number of bases: " + str(baseCount))
    prcntUnkown = (rCount+yCount+otherCount)/baseCount
    print("Percent of unkown bases: "+str(prcntUnkown))







def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument( 
            'files',
            metavar='N',
            type= str,
            nargs='+',
            help='Fasta Files', 
            )
    p.add_argument(
            '-v', '--verbose',
            help='Verbose printing',
            action= 'store_true'
        )
    args = p.parse_args()

    verbose=False


    if args.files == None:
        p.print_help()
        return

    if args.verbose:
        verbose=True


    fastaInfo(args.files, verbose)





if __name__ == '__main__':
    main(sys.argv[1:])

