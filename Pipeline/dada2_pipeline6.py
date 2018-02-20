#!/usr/bin/env python3

##########################################################################################################
############### DADA2 PIPELINE 6a : Convert CSV to FASTA #################################################
##########################################################################################################
####                                                                                                  ####
#### Andreas Novotny, 2018-02                                                                         ####
#### https://github.com/andreasnovotny/DadaSlurm                                                      ####
####                                                                                                  ####
##########################################################################################################
##########################################################################################################
##########################################################################################################

print("Python will now convert CSV sequences to FASTA... ...")

def CSV_to_fasta (input_file, output_file):
	counter = 0
	with open (input_file, 'r') as input_file, open (output_file, 'w') as output_file:
		for row in input_file:
			if counter == 0:
				counter += 1
			else:
				row_list=row.split(',')
				sequence=row_list[1]
				sequence=sequence.replace('\"','')
				x = ''.join(['>', sequence])
				output_file.write(x)
				output_file.write(sequence)
				output_file.write('\n')

if __name__ == '__main__':
	import sys
	filepath=sys.argv[1]
	input_file=''.join([filepath,'/seqs.csv'])
	output_file=''.join([filepath,'/seqs.fa'])

	CSV_to_fasta(input_file, output_file)


##########################################################################################################
##########################################################################################################
