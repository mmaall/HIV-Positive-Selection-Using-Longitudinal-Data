#!/bin/bash/python



from Patient import Patient
import sys
import argparse
import random


class KaKsCalculation: 
	def __init__ (self, patientList=None):

		self.kaksRatios= []
		self.freqTransition= 0
		self.freqTransversion= 0
		if patientList== None:
			self.patientList= []
		else:
			self.patientList= patientList



	def calculateRatio(self):
		for currPatient in PatientList:
			for mutations in currPatient.mutCharDict:
				print(mutations)



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

    print(str(len(patientList)))
    
    random.seed= 1
    currPatientList= []
    for i in range(10):
    	currPatientList.append(patientList[random.randint(0,len(patientList))])

    kaks= KaKsCalculation(currPatientList)




if __name__ == '__main__':
    main(sys.argv[1:])












