#!/usr/bin/env python3
"""
Characterization of mutations between sequences
"""

from MutationCharacterization import FastaReader, mutCharStorage, PatientProfile

class KaKsCalculations:

	def __init__(self, mutDict, patientSize):
		
		self.mutationDict = mutDict
		self.sampleSize = patientSize
		self.sequenceLen = 1143
		
		# all possible comb. of mutations across sample
		self.Ytmut = 0
		self.Yvmut = 0
		self.Stmut = 0
		self.Svmut = 0

		for key in self.mutationDict.keys():
			self.Ytmut += self.mutationDict[key][1]
			self.Yvmut += self.mutationDict[key][2]
			self.Stmut += self.mutationDict[key][3]
			self.Svmut += self.mutationDict[key][4]
		
		# transition/transversion frequencies
		self.t_Freq = (self.Ytmut + self.Stmut) / (self.sequenceLen * self.sampleSize)
		self.v_Freq = (self.Yvmut + self.Svmut) / (2*self.sequenceLen * self.sampleSize)

	def KaKsCalc(self, position, codonInfoList):
		"""
		Performs the actual KaKs calculation.
		"""
		yMutAtCodon = codonInfoList[0]
		sMutAtCodon = codonInfoList[1]
		KaKsval = float(0)
		if self.t_Freq == 0 or self.v_Freq == 0 or sMutAtCodon == 0:
			pass
		else:
			numerator = yMutAtCodon/sMutAtCodon
			bottomDenom = (self.Stmut * self.t_Freq) / (self.Svmut * self.v_Freq)
			topDenom = (self.Ytmut*self.t_Freq) + (self.Yvmut*self.v_Freq)
			fullDenom = topDenom / bottomDenom
			KaKsval = numerator / fullDenom
		return KaKsval
		





import glob, fileinput
from collections import OrderedDict
def main():
	file_list = glob.glob('*.fasta')
	patient_data = []
	myReader = FastaReader()
	patient_list = []
	total_mutation_database = {} # {key: [codon, Yt, Yv, St, Sv]}
	codon_database = {}

	t0 = 0
	tf = 1

	for file in file_list:
		myReader = FastaReader(file)
		for header, sequence in myReader.readFasta():
			patient_data.append([header, sequence])
		mutationStorage = mutCharStorage(patient_data[0][1], patient_data[1][1])
		totalMutations = mutationStorage.findMutations() # Key = Position, value = [codon, Yt, Yv, St, Sv]
		thisPatient = PatientProfile(patient_data[1][0], totalMutations)
		patient_data = [] # resets patient_data
		"""
		UPDATE REQUIRED:
		Just updated the mutation dictionary in MutationCharacterization.py to have this format:
		key: [codon, Yt, Yf, St, Sv]
		Need to make sure that I'm saving the right patient info... must update this loop
		to search for the elements in the format from before.
		"""
		for key in thisPatient.mutations.keys():
			# if the position has not yet been recorded 
			if key not in total_mutation_database.keys():
				num_Yt = 0
				num_Yv = 0
				num_St = 0
				num_Sv = 0
				mutInfo = thisPatient.mutations[key]
				for item in range(len(mutInfo)):
					if mutInfo[1] == 1:
						num_Yt += 1
					if mutInfo[2] == 1:
						num_Yv += 1
					if mutInfo[3] == 1:
						num_St += 1
					if mutInfo[4] == 0:
						num_Sv = 0

				total_mutation_database.update({key:[mutInfo[0], num_Yt, num_Yv, num_St, num_Sv]})
			else:
				mutInfo = thisPatient.mutations[key]
				for item in range(len(mutInfo)):
					if mutInfo[1] == 1:
						total_mutation_database[key][1] += 1
					if mutInfo[2] == 1:
						total_mutation_database[key][2] += 1
					if mutInfo[3] == 1: 
						total_mutation_database[key][3] += 1
					if mutInfo[4] == 1: 
						total_mutation_database[key][4] += 1
	
	# makes dictionary for individual codons
	for key in total_mutation_database.keys():
		thisCodon = total_mutation_database[key][0]
		if thisCodon not in codon_database:

			codon_database.update({thisCodon: [(total_mutation_database[key][1]+total_mutation_database[key][2]), (total_mutation_database[key][3]+total_mutation_database[key][4])]}) # syn, nonsyn
		else:
			codon_database[thisCodon][0] += (total_mutation_database[key][1]+total_mutation_database[key][2]) # adds syn
			codon_database[thisCodon][1] += (total_mutation_database[key][3]+total_mutation_database[key][4])

	allMutKaKs = KaKsCalculations(total_mutation_database, len(file_list))
	total_mutation_database = OrderedDict(sorted(total_mutation_database.items(), key=lambda kv: kv[1]))
	
	for key in total_mutation_database.keys():
		originalMutInfo = total_mutation_database[key]
		codon = total_mutation_database[key][0]
		codonInfo = codon_database[codon]
		KaKsValue = allMutKaKs.KaKsCalc(key, codonInfo)
		originalMutInfo.append(KaKsValue)
		print("Mutation "+str(key)+" @ codon "+str(total_mutation_database[key][0])+": "+str(KaKsValue))




if __name__ == '__main__':
	main()

	

