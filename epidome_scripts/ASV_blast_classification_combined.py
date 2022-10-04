#!/usr/bin/env python3

import os
import sys
import argparse

def set_up_argparse():
	parser = argparse.ArgumentParser(description="""
	**
	Classify dada ASVs using blastn similarity.
	Requires dada output for both primer pairs and reference fasta files for both.
	""")

	parser.add_argument('--p1_file', '-p1', help='Dada ASV output for primer 1 (g216)')
	parser.add_argument('--p2_file', '-p2', help='Dada ASV output for primer 2 (yycH)')
	parser.add_argument('--p1_ref', '-p1_ref', help='Reference fasta for primer 1 (g216)')
	parser.add_argument('--p2_ref', '-p2_ref', help='Reference fasta for primer 2 (yycH)')
	parser.add_argument('--pident_req', '-pident', help='Minimum percent identity to classify ASV',default='99.5')
	opts = parser.parse_args()
	
	return opts

def check_if_flag_is_provided(flag,name,option):
	if flag is None:
		print("Checking if", name, "is provided...No!")
		print("You have not provided the", name,"!", "use the",option,"option to provide the", name)
		sys.exit(1)
	else:
		print("Checking if", name, "is provided...Yes")


def setup_fasta(dada_csv_file,fasta_output_file):
	csv_lines = []
	firstlineflag = 1
	fasta_print = ''
	with open(dada_csv_file) as f:
		for line in f:
			line = line.rstrip('\n')
			if firstlineflag==1:
				firstlineflag = 0
				header = line
			else:
				splitline = line.split(";")
				csv_lines.append(splitline)
				seq_number = splitline[0][1:-1]
				seq = splitline[1][1:-1]
				asv_name = splitline[2]
				fasta_print += '>'+seq_number+'\n'+seq+'\n'
	o = open(fasta_output_file,'w')
	o.write(fasta_print)
	o.close()
	return(header, csv_lines)

def load_dada_output(dada_csv_file):
	firstlineflag = 1
	counts = []
	IDs = []
	with open(dada_csv_file) as f:
		for line in f:
			line = line.rstrip('\n')
			if firstlineflag == 1:
				firstlineflag = 0
				seqs = line.split(';')[1:]
			else:
				splitline = line.split(';')
				counts.append(splitline[1:])
				IDs.append(splitline[0])
	return(seqs,IDs,counts)

def print_fasta(seqs,ASV_fasta):
	printline = ''
	n=0
	for seq in seqs:
		n+=1
		printline += '>'+str(n)+'\n'+seq[1:-1]+'\n'
	o = open(ASV_fasta,'w')
	o.write(printline)
	o.close()

def run_blast(ASV_fasta,reference_fasta,blast_output_file,pident_req):
	cmd = 'blastn -query '+ASV_fasta+' -subject '+reference_fasta+' -out '+blast_output_file+' -outfmt 6 -perc_identity '+pident_req
	os.system(cmd)

def parse_blast(blast_output_file,pident_req):
	ASV_to_ref = {}
	bitscores = {}
	pident_req = float(pident_req)
	with open(blast_output_file) as f:
		for line in f:
			line = line.rstrip('\n').split('\t')
			pident = float(line[2])
			ASV_number = line[0]
			ref_number = line[1]
			bitscore = float(line[11])
			if ASV_number in ASV_to_ref:
				if bitscores[ASV_number]<bitscore:
					ASV_to_ref[ASV_number] = [ref_number]
					bitscores[ASV_number] = bitscore
				elif bitscores[ASV_number]==bitscore:
					ASV_to_ref[ASV_number].append(ref_number)
			else:
				ASV_to_ref[ASV_number] = [ref_number]
				bitscores[ASV_number] = bitscore
	return(ASV_to_ref)

def rename_ASVs(ASV_to_ref_dict,csv_header,csv_lines,dada_csv_renamed_file):
	printline = csv_header+'\n'
	for line in csv_lines:
		ASV_number = line[0][1:-1]
		if ASV_number in ASV_to_ref_dict:
			new_ASV_number = "\"seq"+ASV_to_ref_dict[ASV_number]+"\""
			old_ASV_number = line[2]
			if new_ASV_number != old_ASV_number:
				print(old_ASV_number+' renamed to '+new_ASV_number, file=sys.stderr)
				line[2] = new_ASV_number
		printline += ';'.join(line)+'\n'
	o = open(dada_csv_renamed_file,'w')
	o.write(printline)
	o.close()

def print_ASV_count_table(seqs, IDs, counts, ASV_to_ref,new_dada_output_file):
	printline = "\"ASV\";\"Seq_number\";"+";".join(IDs)+"\n"
	for n in range(len(seqs)):
		ASV = seqs[n][1:-1]
		count_list = []
		for m in range(len(counts)):
			count_list.append(counts[m][n])
		number = str(n+1)
		if number in ASV_to_ref:
			seq_number = ','.join(ASV_to_ref[number])
			if seq_number == "NA":
				printline += "\""+str(number)+"\";\""+ASV+"\";"+seq_number+";"+";".join(count_list)+"\n"
			else:
				printline += "\""+str(number)+"\";\""+ASV+"\";\"seq"+seq_number+"\";"+";".join(count_list)+"\n"
		else:
			seq_number = "NA"
			printline += "\""+str(number)+"\";\""+ASV+"\";"+seq_number+";"+";".join(count_list)+"\n"
	o = open(new_dada_output_file,'w')
	o.write(printline)
	o.close()


def main(opts):
	check_if_flag_is_provided(opts.p1_file,"primer 1 (g216)","-p1")
	check_if_flag_is_provided(opts.p2_file,"primer 2 (yycH)","-p2")
	check_if_flag_is_provided(opts.p1_ref,"primer 1 reference fasta (g216)","-p1_ref")
	check_if_flag_is_provided(opts.p2_ref,"primer 2 reference fasta (yycH)","-p2_ref")
	#outdir = os.path.split(opts.p1_file)[0]
	p1_fasta = opts.p1_file+'.ASV_seqs.fasta'
	p2_fasta = opts.p2_file+'.ASV_seqs.fasta'
	p1_blast = opts.p1_file+'.ASV_blast.txt'
	p2_blast = opts.p2_file+'.ASV_blast.txt'
	p1_out = opts.p1_file+'.classified.csv'
	p2_out = opts.p2_file+'.classified.csv'
	print("Classifying primer 1 sequences (g216)")
	p1_seqs,p1_IDs,p1_counts = load_dada_output(opts.p1_file)
	print_fasta(p1_seqs,p1_fasta)
	run_blast(p1_fasta,opts.p1_ref,p1_blast,opts.pident_req)
	ASV_to_ref_p1 = parse_blast(p1_blast,opts.pident_req)
	print_ASV_count_table(p1_seqs,p1_IDs,p1_counts,ASV_to_ref_p1,p1_out)
	print("Classifying primer 2 sequences (yycH)")
	p2_seqs,p2_IDs,p2_counts = load_dada_output(opts.p2_file)
	print_fasta(p2_seqs,p2_fasta)
	run_blast(p2_fasta,opts.p2_ref,p2_blast,opts.pident_req)
	ASV_to_ref_p2 = parse_blast(p2_blast,opts.pident_req)
	print_ASV_count_table(p2_seqs,p2_IDs,p2_counts,ASV_to_ref_p2,p2_out)
	print("Classified primer 1 sequences printed to "+p1_out)
	print("Classified primer 2 sequences printed to "+p2_out)
	return;


if __name__ == '__main__':
	opts = set_up_argparse()
	main(opts)