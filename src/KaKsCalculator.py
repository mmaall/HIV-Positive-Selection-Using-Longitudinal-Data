#!/bin/bash/python



from Patient import Patient
from Patient import transTranv
from Patient import mutationType
from Patient import baseType
import sys
import argparse
import random


class KaKsCalculation: 
	def __init__ (self, patientList=None):
		self.kaksRatioByBase= []
		self.kaksRatioByCodon= []
		self.freqTransition= 0
		self.freqTransversion= 0
		self.standardSeqLength= 1680
		if patientList== None:
			self.patientList= []
		else:
			self.patientList= patientList

	def kaksHelper(synCount, nonSynCount, transFreq, transvFreq, possMutationList):
		numerator= nonSynCount/ synCount
		topDenominator= (transFreq* possMutationList[1]) + (transvFreq* possMutationList[3])
		botDenominator=	(transFreq* possMutationList[0]) + (transvFreq* possMutationList[2])
		return numerator / (topDenominator/botDenominator)
	def calculateRatio(self):
		#We're not going to use the 1s place for these
		totalNucleotides=0
		numPatients= len(self.patientList)
		seqLength= self.standardSeqLength +1
		#List containing the number of synonymous mutations at a position
		synCount= [0] * seqLength
		#List containing the number of nonsynonymous mutations at a given position
		nonSynCount= [0] * seqLength
		#List containing the counts of possible mutations at a specific positions	
		#Every element is a list of four elements 
		#[synonymous transitions, nonsynonymous transitions, synonymous transversions, nonsynonymous transversions]	
		possMutations= [[0,0,0,0]]* seqLength
		transCount= 0
		transvCount= 0

		for currPatient in self.patientList:
			totalNucleotides= currPatient.effectiveLength
			for key in currPatient.mutCharDict:
				# Key = Position, value = [t0 base, tf base, transition/transversion, mutType]
				mutationInfo = currPatient.mutCharDict[key]
				if mutationInfo[2] == transTranv.transition:
					#transCounts[key]+=1
					transCount+=1
				else: 
					#transvCounts[key]+=1
					transvCount+=1

				if mutationInfo[3] == mutationType.syn:
					synCount[key]+=1
				else:
					nonSynCount[key]+=1

			currIndex= 0
			for mutationTuple in currPatient.possibleMutations:
				if mutationTuple != None:
					synTrans, nonSynTrans, synTransv, nonSynTransv = mutationTuple
					#print(mutationTuple)
					#print(str(possMutations[currIndex][0]))
					possMutations[currIndex][0]+=synTrans
					possMutations[currIndex][1]+=nonSynTrans
					possMutations[currIndex][2]+=synTransv
					possMutations[currIndex][3]+=nonSynTransv
				currIndex+=1
		#Calculate Transition/transversion frequency
		transFreq= transCount/(totalNucleotides*numPatients)
		transvFreq= transvCount/(totalNucleotides*numPatients*2)
		for numSyn, numNonSyn, mutList in zip(synCount[1:], nonSynCount[1:], possMutations[1:]):
			ratio=0
			if  numSyn == 0:
				#numSyn=1
				ratio = 0
				self.kaksRatioByBase.append(ratio)
				continue
			numerator= numNonSyn/ numSyn
			topDenominator= (transFreq* mutList[1]) + (transvFreq* mutList[3])
			botDenominator=	(transFreq* mutList[0]) + (transvFreq* mutList[2])
			ratio= numerator / (topDenominator/botDenominator)
			self.kaksRatioByBase.append(ratio)
			#self.kaksRatioByBase.append(self.kaksHelper(numSyn, numNonSyn,transFreq,transvFreq ))




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

    patientList= []

    for fileName in args.files:
    	newPatient= Patient()
    	newPatient.inputFile(fileName)
    	patientList.append(newPatient)

    #print(str(len(patientList)))
    
    print("Patients Inputed")
    random.seed= 1


    #Holds the averages for the Kaks ratios
    kaksAvg= [0] * (1680+1)
    #pvalues
    pValues= [0] * (1680+1)
    #Holds the averages for the 
    numReplicates= 10000 # Number of bootstrap replicates
    numPatients= len(patientList)
    selectionPct= .2
    patientsPerBootstrap= int(selectionPct*numPatients)

    for i in range(numReplicates):
    	print("Replicate: " + str(i))
    	sampleList= random.sample(range(0, numPatients), patientsPerBootstrap)
    	selectedPatientList =[]
    	for patientIndex in sampleList:
    		selectedPatientList.append(patientList[patientIndex])
    	currKaKs= KaKsCalculation(selectedPatientList) 
    	currKaKs.calculateRatio()
    	for position in range(0,len(currKaKs.kaksRatioByBase)):
    		kaksAvg[position]+= currKaKs.kaksRatioByBase[position]
    for i in range(len(kaksAvg)):
    	kaksAvg[i]= kaksAvg[i]/numReplicates

    count= 0
    for ratio in kaksAvg:
    	if(ratio>= 1.0):
    		codon= int((count-1)/3)+1
    		print("Position: "+str(count)+"\tKaKs Ratio: "+str(ratio)+"\tCodon: " +str(codon))
    	count +=1
  



if __name__ == '__main__':
    main(sys.argv[1:])












