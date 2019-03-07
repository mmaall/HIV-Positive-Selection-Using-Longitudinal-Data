#!/bin/bash/python



from Patient import Patient
import sys
import argparse



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
    	print("Inputing patient "+fileName)
    	newPatient= Patient()
    	newPatient.inputFile(fileName)
    	patientList.append(newPatient)

    print(str(len(patientList)))
    for patient in patientList:
    	print (patient)




if __name__ == '__main__':
    main(sys.argv[1:])












