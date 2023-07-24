#!/usr/bin/env python

"""
This script performs a quality control of orthologous protein-coding sequences, which must be arrenged in the correct frame. 
Filter out sequence based on length variation.
Remove distantly related sequences based on a predefined cutoff. 
Remove poorly informative sequences based on a predifined threshold of the proportion of gaps and missing data.
"""

import argparse,os,sys
import logging
import csv
import re
import subprocess
import numpy
import shutil
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO



#parser=MyParser()
parser = argparse.ArgumentParser(description="Filter out distantly related samples from fasta files of DNA sequences of cluster of orthologs.")
parser.add_argument('--fasta_file', help='Path of fasta file.')
#parser.add_argument('--p', help='Number of processors')
parser.add_argument('--log_file', help='Path to the log file')
parser.add_argument('--out_dir', help='Path to the output dir')
parser.add_argument('--distance_threshold', help='Maximum average of uncorrected pairwise distance per sample. Possible values from 0 to 100.')
parser.add_argument('--proportion_Ns_gaps', help='Maximum fraction of allowed missing data and gaps. For example, 0.1 means that sequences with more than 10%% of missing data and gaps will be deleted.')
parser.add_argument('--length_filter', help='Cutoff of length variation. For example, 10 means that sequences with a length greater than 110%% and lower than 90%% of the average length will be deleted.')
parser.add_argument('--cleanup', help='If it sets to yes, all intermediate files will be deteled.')


if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)

args = parser.parse_args()

if args.fasta_file:
	fasta_file = args.fasta_file

#if args.p:
#	num_cores = args.p

if args.log_file:
	log_file = args.log_file

if args.out_dir:
	out_dir = args.out_dir
	if out_dir[-1] != "/":
		out_dir += "/"

if args.distance_threshold:
	distance_threshold = float(args.distance_threshold)

if args.proportion_Ns_gaps:
	proportion_Ns_gaps = float(args.proportion_Ns_gaps)

if args.cleanup:
	cleanup = args.cleanup

if args.length_filter:
	length_filter = float(args.length_filter)

"""
Functions
"""



def create_dir(dir_path):
	"""Create a new directory
	"""
	if not os.path.exists(dir_path):
		os.makedirs(dir_path)
	return dir_path

def remove_by_stop_m3_size_trim_Ns(path_fasta,length_filter):
	"""Remove Ns at the edges and gaps in any position.
	Remove sequences not multiple of three, sequences that contain premature stop codons and sequences based on length variation.
	"""
	#Create the lists
	list_remove_multiple_3 = []
	list_remove_stop_codon = []
	list_remove_by_length = []
	lengths = []
	ids = []

	#Create a dictionary of records using Biopython
	records = SeqIO.to_dict(SeqIO.parse(path_fasta, "fasta"))
	
	#Number o sequences
	original_num_seq = len(records)
	
	#Iterate in the alignment.
	for record_index in records:
		seq = str(records[record_index].seq)
		
		#Remove Ns at the edges and all gaps using the function trim_Ns_gaps().
		seq = trim_Ns_gaps(seq)
		records[record_index].seq = Seq(seq)

		seq_length = len(seq)
		#Get sequences not multiple of 3.
		if seq_length % 3 != 0:
			list_remove_multiple_3.append(records[record_index].id)

		else:
			#Get sequences with a premature stop codon.
			for pos_codon in range(0,seq_length -3,3):
				codon = seq[pos_codon]+seq[pos_codon+1]+seq[pos_codon+2]
				if  codon == "TAA" or codon == "TAG" or codon == "TGA":
					list_remove_stop_codon.append(records[record_index].id)
					break

	#Create a list containing all sequences to be removed.
	list_all_to_eliminate = list_remove_multiple_3 + list_remove_stop_codon
	
	#Remove sequences
	for id in list_all_to_eliminate:
		if id in records:
			records.pop(id)
	
	#Convert dictionary of records to a list of records:
	list_records = []
	for key in records:
		list_records.append(records[key])

	#Compare lengths
	for record in list_records:
		#Get the length of the sequences
		lengths.append(len(record.seq))
		ids.append(record.id)


	#Remove sequences based on a cutoff value of length variation
	ids_lengths = [ids,lengths]

	list_remove_by_length_index = filter_by_length(ids_lengths,list_remove_by_length,length_filter,list_records)

	#Filter list of records by length
	list_remove_by_length_index = sorted(list_remove_by_length_index, reverse=True)
	for idx in list_remove_by_length_index:
		if idx < len(list_records):
			list_records.pop(idx)
	
	return list_remove_multiple_3,list_remove_stop_codon,list_remove_by_length,original_num_seq,list_records

def trim_Ns_gaps(sequence):
	"""Trim missing data from the ends of a sequence and delete all gaps (-). The sequence must be a string.
	"""

	#Convert to uppercase
	new_seq = sequence.upper()
	#Trim the 5 prime end
	five_prime_trim = re.sub("^[N]*","",new_seq)
	#Trim the 3 prime end
	trimmed_seq = re.sub("[N]*$","",five_prime_trim)
	#Delete all gaps
	clean_seq = trimmed_seq.replace("-", "")
	return clean_seq


def filter_by_length(ids_lengths,list_remove_by_length,length_filter,list_records):
	"""Generate a list of record IDs to be removed based on the length variation cutoff.
	"""
	list_remove_by_length_index = []
	eliminate_tmp = ["tmp"]
	while len(eliminate_tmp) > 0:
		eliminate_tmp = []
		lengths = ids_lengths[1]
		Maximum_length, Minimum_length = calculate_max_min_length(lengths,length_filter)
		max_val,max_index,min_val,min_index = max_min(lengths)

		diff_max = 0
		diff_min = 0
		if max_val > Maximum_length:
			diff_max = max_val - Maximum_length

		if min_val < Minimum_length:
			diff_min = Minimum_length - min_val

		if diff_max == 0 and diff_min == 0:
			break
		
		else:
			if diff_max >= diff_min:
				eliminate_tmp.append(ids_lengths[0][max_index])
				list_remove_by_length.append(ids_lengths[0][max_index])
				ids_lengths[0].pop(max_index)
				ids_lengths[1].pop(max_index)
			else:
				eliminate_tmp.append(ids_lengths[0][min_index])
				list_remove_by_length.append(ids_lengths[0][min_index])
				ids_lengths[0].pop(min_index)
				ids_lengths[1].pop(min_index)

	for record_index in range(len(list_records)):
		for id_remove in list_remove_by_length:
			if list_records[record_index].id == id_remove:
				list_remove_by_length_index.append(record_index)
				break

	return list_remove_by_length_index

 
def max_min(list):
	"""Retrieve the maximum and minimum values of a list and their indexes.
	"""
	max_val = max(list)
	max_index = list.index(max_val)
	min_val = min(list)
	min_index = list.index(min_val)
	return max_val,max_index,min_val,min_index


def calculate_max_min_length(lengths,length_filter):
	"""Calulate the the maximum and minimum value
	"""
	number_seqs = len(lengths)
	average = round(float(sum(lengths))/number_seqs,2)

	Maximum_length = average*(100+length_filter)/100
	Minimum_length = average*(100-length_filter)/100
	return Maximum_length, Minimum_length


def save_results(records,output_path):
	"""Save a Biopython list of records into a file in fasta format.
	"""
	with open(output_path, "w") as output_handle:
		for id in range(len(records)):
			output_handle.write(">" + records[id].id + "\n" + str(records[id].seq) + "\n") 
	return output_path


def filter_distant_seqs(list_records_dna,list_ids,distance_matrix_array,distance_threshold):
	"""#Filter out dna sequences based on a distance matrix
	"""
	list_removed_ids = []
	list_keep_ids = []
	backup_list_records_dna = list_records_dna[:]
	#Iterate in the distance matrix
	for index in range(len(distance_matrix_array)):
		sample_dist_list = distance_matrix_array[index]
		#The number of elements, we need to eliminate the zero (diagonal value)
		number_elements = len(sample_dist_list) -1
		average = round(float(sum(sample_dist_list))/number_elements,2)
		if average > distance_threshold:
			for record in backup_list_records_dna:
				if record.id == list_ids[index]:
					list_removed_ids.append(record.id)
					break
		else:
			for record in backup_list_records_dna:
				if record.id == list_ids[index]:
					list_keep_ids.append(record.id)
					break
	
	return list_keep_ids,list_removed_ids

def prepare_colnumbering_trimal(stdout):
	"""Parse the standard output from the option -colnumbering of trimal and convert it into a list of integers.
	"""
	column_list=str(stdout).split(", ")
	column_list[0]=column_list[0][4:]
	column_list[-1]=column_list[-1][:-5]
	
	for index in range(len(column_list)):
		column_list[index] = int(column_list[index])

	return column_list

def get_sites_dna_from_trimal_aa_output(translatorX_nuc_out,column_info_list):
	"""Edit the alignment in codons to. This function saves the results in a fasta file and returns a biopython list of records.
	"""

	output_trimal_dna = translatorX_nuc_out.split(".fasta")[0] + "_strictplus.fasta"
	aln_dna = AlignIO.read(open(translatorX_nuc_out), 'fasta')

	#Use the column info from trimal output of aa alignemnt to create a DNA aligment (codon based)
	for record_index in range(len(aln_dna)):
		sequence = str(aln_dna[record_index].seq)
		new_sequence = ""
		for index_aa in column_info_list:
			index_dna = index_aa*3
			new_sequence = new_sequence + sequence[index_dna:index_dna+3]
		aln_dna[record_index].seq = Seq(new_sequence)
	
	#Save the results in a fasta file
	save_results(aln_dna,output_trimal_dna)
			
	return output_trimal_dna,aln_dna

def align_translatorx(fasta_path,output_prefix,output_cluster_dir):
	"""Run the translator X software to align the sequences using muscle algorithm. 
	This function returns the path of the alignments in aa and dna. Temporary files will be deleted.
	"""
	translatorX_nuc_out = output_cluster_dir + output_prefix + ".nt_ali.fasta"
	translatorX_aa_out = output_cluster_dir + output_prefix + ".aa_ali.fasta"

	cmd = ["translatorx", "-i",fasta_path,"-o",output_prefix]
	p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=output_cluster_dir)
	out = p.communicate()

	#Path of secondary outputs
	html_coloured = output_cluster_dir + output_prefix + ".aa_based_codon_coloured.html"
	html = output_cluster_dir + output_prefix + ".html"
	aaseqs = output_cluster_dir + output_prefix + ".aaseqs"
	aaseqs_fasta = output_cluster_dir + output_prefix + ".aaseqs.fasta"
	muscle_log = output_cluster_dir + output_prefix + ".muscle.log"
	nt12_ali = output_cluster_dir + output_prefix + ".nt12_ali.fasta"
	nt1_ali = output_cluster_dir + output_prefix + ".nt1_ali.fasta"
	nt2_ali = output_cluster_dir + output_prefix + ".nt2_ali.fasta"
	nt3_ali = output_cluster_dir + output_prefix + ".nt3_ali.fasta"
	
	#Remove intermediate files
	try:
		os.remove(html_coloured)
		os.remove(html)
		os.remove(aaseqs)
		os.remove(aaseqs_fasta)
		os.remove(muscle_log)
		os.remove(nt12_ali)
		os.remove(nt1_ali)
		os.remove(nt2_ali)
		os.remove(nt3_ali)

	except: pass # no need to worry about extra intermediate files

	return translatorX_nuc_out,translatorX_aa_out


def trimal_stricplus(translatorX_aa_out,output_cluster_dir):
	"""Run trimal to filter the protein alignment.
	"""

	output_trimal = translatorX_aa_out.split(".fasta")[0] + "_strictplus.fasta"
	cmd = ["trimal","-strictplus","-colnumbering","-in",translatorX_aa_out,"-out",output_trimal]
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd= output_cluster_dir)
	column_out = p.stdout.read()
	out = p.communicate()
	column_list = prepare_colnumbering_trimal(column_out)	
	return column_list


def prepare_distmat_remove_Ns(aln_dna):
	"""Remove all sites with missing data (Ns) from an alignment. 
	The output of this function will only be used as the input of the Distmat program of EMBOSS toolkit (Rice et al. 2000).
	"""
	#Creates a temporary file
	output_trimal_dna_temp = aln_dna.split(".fasta")[0] + "_temp.fasta"
	aln_dna = AlignIO.read(open(aln_dna), 'fasta')
	
	#Get the indexes for all Ns
	index_final_seq=[]
	index_Ns=[]
	for record in aln_dna:
		seq = str(record.seq)
		for position in range(len(seq)):
			if seq[position] == "N":
				if position not in index_Ns:
					index_Ns.append(position)

	#Get the sequences wihtout Ns
	for record_index in range(len(aln_dna)):
		seq = str(aln_dna[record_index].seq)
		new_seq = ""
		for index in range(len(seq)):
			if index not in index_Ns:
				new_seq = new_seq + seq[index]
		aln_dna[record_index].seq = Seq(new_seq)

	#Save the results in a fasta file
	save_results(aln_dna,output_trimal_dna_temp)
	return output_trimal_dna_temp

def retrieve_fasta_dict(id_list,input_fasta,output_fasta):
	"""Retrieve a set of sequences in fasta format file based on a list of IDs and a multifasta file, which contains the record of these IDs.
	"""
	record_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
	with open(output_fasta, "w") as output_handle:
		for ID_seq_record in id_list:
			if ID_seq_record in record_dict:
				output_handle.write(record_dict[ID_seq_record].format("fasta"))
	return output_fasta


def distmat(output_trimal_dna,output_cluster_dir):
	"""Convert the distmat result to a distance matrix in phylip format
	"""
	#Create the path to trimal dna dmat
	output_trimal_dna_dmat = output_trimal_dna.split(".fasta")[0] + ".dist"

	#Run distmat
	cmd = ["distmat","-nucmethod","0","-sequence",output_trimal_dna,"-outfile",output_trimal_dna_dmat] 
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, cwd= output_cluster_dir)
	out = p.communicate()

	return output_trimal_dna_dmat


def distmat_to_np_array(output_trimal_dna_dmat):
	"""Parse distmat output to a numpy array
	"""

	#Open the dist mat file
	with open(output_trimal_dna_dmat,"r") as dist_handle:
		lines_dist_matrix = dist_handle.readlines()

	#Create the list to put the results
	distance_matrix_final=[]
	list_ids = []

	#Populate the list of distance matrix and list of ids
	for index_dist in range(8,len(lines_dist_matrix)):
		list_dist = lines_dist_matrix[index_dist].split("\t")[1:-2]
		list_ids.append(lines_dist_matrix[index_dist].split("\t")[-1].strip().split(" ")[0])
		for index_pair_dist in range(len(list_dist)):
			if list_dist[index_pair_dist] == '':
				list_dist[index_pair_dist] = float(0.00)
			else:
				list_dist[index_pair_dist] = float(list_dist[index_pair_dist])
		distance_matrix_final.append(list_dist)

	#Convert into a numpy array
	distance_matrix_array = numpy.array(distance_matrix_final)
	#Transpose the uuper diagonal matrix and convert it into a complete matrix
	distance_matrix_array = distance_matrix_array.T + distance_matrix_array
	return list_ids,distance_matrix_array



def wraper_translatorx_trimal(path_new_fasta,translatorX_prefix,path_output_directory):
	"""Align and filter the alignment based on aa.
	"""
	#Align protein sequences using translatorX
	translatorX_nuc_out,translatorX_aa_out = align_translatorx(path_new_fasta,translatorX_prefix,path_output_directory)
	#Filter the aligment
	column_info = trimal_stricplus(translatorX_aa_out,path_output_directory)
	#Get the nucleotide aligment based on aa alignment filtered by trimal
	output_trimal_dna,list_records_dna = get_sites_dna_from_trimal_aa_output(translatorX_nuc_out,column_info)

	return output_trimal_dna,list_records_dna

#Save info from a list in the log file
def list_to_log_file(list, reason, clusterID):
	"""Save info from a list in the log file
	"""
	if len(list)>0:
		for sample in list:
			logging.info(clusterID + ": " + sample + " has been removed due to " + reason)
		return clusterID
	else:
		return clusterID


def filter_by_Ns_gaps(fasta_path,maximum_proportion_Ns_gaps):
	"""The function removes sequences with less than an arbitrary proportion of missing data (Ns) and gaps.
	"""
	output_Ns_gaps = fasta_path.split(".fasta")[0] + "_Ns_gaps.fasta"
	record_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
	#number_seq = len(record_dict) 
	records = []
	list_remove_Ns_gaps = []
	list_record_ids = []
	for key in record_dict:
		record = record_dict[key]
		sequence = str(record.seq).upper()
		length_alignment = len(sequence)
		count_dict = Counter(sequence)
		count_Ns_gaps = count_dict["N"] + count_dict["-"]
		Proportion_Ns_gaps = round(count_Ns_gaps/length_alignment,2)
		
		#Compare with an arbitrary value
		if Proportion_Ns_gaps <= maximum_proportion_Ns_gaps:
			records.append(record)
			list_record_ids.append(record.id)
		else:
			list_remove_Ns_gaps.append(record.id)
	
	#Save the results
	if len(list_remove_Ns_gaps) != 0:
		save_results(records,output_Ns_gaps)
	return list_record_ids,list_remove_Ns_gaps,output_Ns_gaps

def cleaning_intermediate_files(output_dir,ClusterID):
	"""This functions removes intermediate files
	"""
	try:
		os.remove(output_dir + ClusterID + ".aa_ali.fasta")
		os.remove(output_dir + ClusterID + ".aa_ali_strictplus.fasta")
		os.remove(output_dir + ClusterID + ".nt_ali.fasta")
		os.remove(output_dir + ClusterID + ".nt_ali_strictplus.fasta")
		os.remove(output_dir + ClusterID + ".nt_ali_strictplus_temp.dist")
		os.remove(output_dir + ClusterID + ".nt_ali_strictplus_temp.fasta")
		os.remove(output_dir + ClusterID + "_filter_dist.fas")
		os.remove(output_dir + ClusterID + "_filter_dm.aa_ali.fasta")
		os.remove(output_dir + ClusterID + "_filter_dm.aa_ali_strictplus.fasta")
		os.remove(output_dir + ClusterID + "_filter_dm.nt_ali.fasta")
		os.remove(output_dir + ClusterID + "_filter_dm.nt_ali_strictplus.fasta")
		os.remove(output_dir + ClusterID + "_filter_dm_Ns_gaps.aa_ali.fasta")
		os.remove(output_dir + ClusterID + "_filter_dm_Ns_gaps.aa_ali_strictplus.fasta")
		os.remove(output_dir + ClusterID + "_filter_dm_Ns_gaps.fas")
		os.remove(output_dir + ClusterID + "_filter_dm_Ns_gaps.nt_ali.fasta")
		os.remove(output_dir + ClusterID + "_filter_dm.nt_ali_strictplus_Ns_gaps.fasta")
	
	except: pass # no need to worry about extra intermediate files

	return ClusterID

def number_seq_fasta(input_fasta):
	"""This function returns the number of sequences of a fasta file
	"""
	record_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
	return len(record_dict)


def cleaning(fasta_file,out_dir,distance_threshold,proportion_Ns_gaps,length_filter,cleanup):
	"""Whrapper function that summarizes all the steps.
	"""
	if os.path.isfile(fasta_file):
	
		#Get the cluster ID and paths
		path_fasta_file = fasta_file
		fasta = os.path.basename(path_fasta_file)
		ClusterID = fasta.split(".")[0]
		path_output_directory = out_dir + ClusterID + "/"
		path_new_fasta = path_output_directory + fasta
		
		#Make the new directory using the clusterID
		create_dir(path_output_directory)
		
		#Filter sequences
		list_remove_multiple_3,list_remove_stop_codon,list_remove_by_length,original_num_seq,records = remove_by_stop_m3_size_trim_Ns(path_fasta_file,length_filter)
		
		number_seq_removed = len(list_remove_multiple_3) + len(list_remove_stop_codon) + len(list_remove_by_length)
		assert number_seq_removed + len(records) == original_num_seq, "WARNING: Something is wrong in the function remove_by_stop_m3_size_trim_Ns."

		#Print the info in the log file
		list_to_log_file(list_remove_multiple_3, "the sequence is not multiple of 3", ClusterID)
		list_to_log_file(list_remove_stop_codon, "the sequence has stop codons", ClusterID)
		list_to_log_file(list_remove_by_length, "the sequence has been removed based on its length.", ClusterID)

		#test if the file is not empty
		if len(records) != 0:
			#Save the results into a fasta file
			save_results(records,path_new_fasta)
		else:
			logging.info(clusterID + ": All sequences have been removed due to they are not multiple of 3 or unexpected stop codons")
			return 
	
		#Align and filter based on trimal
		translatorX_prefix = ClusterID
		output_trimal_dna,list_records_dna = wraper_translatorx_trimal(path_new_fasta,translatorX_prefix,path_output_directory)		
		
		#Remove Ns to run distmat
		output_trimal_dna_temp = prepare_distmat_remove_Ns(output_trimal_dna)
		
		#Calculate pairwise distance using distmat (Emboss)
		output_trimal_dna_dmat = distmat(output_trimal_dna_temp,path_output_directory)
		
		#Parse the distmat results
		list_ids,distance_matrix_array = distmat_to_np_array(output_trimal_dna_dmat)
		#Filtering based on distance matrix
		list_keep_ids,list_removed_ids = filter_distant_seqs(list_records_dna,list_ids,distance_matrix_array,distance_threshold)
		
		assert len(list_keep_ids) + len(list_removed_ids) == len(list_records_dna), "WARNING: Something is wrong in the functions for filtering using the distance matrix."

		#Print a log message with the removed sequences
		log_message = "the sequence has more than " + str(distance_threshold) + " of average genetic distance."
		list_to_log_file(list_removed_ids,log_message,ClusterID)
		
		if len(list_keep_ids) != 0:
			#Get the filtered sequences
			fasta_filtered_by_dist_matrix = path_output_directory + ClusterID + "_filter_dist.fas"
			retrieve_fasta_dict(list_keep_ids,path_new_fasta,fasta_filtered_by_dist_matrix)
		else:
			logging.info(ClusterID + ": All sequences have been removed due to they are not multiple of 3, unexpected stop codons or has more than " + str(distance_threshold) + " of average genetic distance.")
			return

		#Second round filter by missing data (gap and Ns)
		#Do the alignment using translatorX and filter the alignment using trimal
		translatorX_prefix_filtered_dist = ClusterID + "_filter_dm"
		output_trimal_dna,list_records_dna = wraper_translatorx_trimal(fasta_filtered_by_dist_matrix,translatorX_prefix_filtered_dist,path_output_directory)

		#Filter using a proportion of Ns and gaps
		list_record_ids,list_remove_Ns_gaps,output_Ns_gaps = filter_by_Ns_gaps(output_trimal_dna,proportion_Ns_gaps)
		
		#Get the number of sequences
		number_seq = number_seq_fasta(output_trimal_dna)
		assert len(list_record_ids) + len(list_remove_Ns_gaps) == number_seq, "WARNING: Something is wrong in the theat filters Ns and gaps."

		#Print a log message with the removed sequences
		log_message = "the sequence has more than " + str(proportion_Ns_gaps) + " proportion of Ns and gaps."
		list_to_log_file(list_remove_Ns_gaps,log_message,ClusterID)
		
		#test if the file is not empty
		if len(list_record_ids) != 0:
			#Realign if necessary
			if len(list_remove_Ns_gaps) == 0:
				#Save the final filtered fasta file
				final_filtered_fasta = path_output_directory + ClusterID + "_filtered.fas"
				shutil.move(output_trimal_dna,final_filtered_fasta)
				#Remove intermediate files
				if cleanup == "yes":
					cleaning_intermediate_files(path_output_directory,ClusterID)
				
				logging.info(ClusterID + ": Has been successefully filtered.")
				return
			else:
				#Get the filtered sequences
				fasta_filtered_by_dist_matrix_Ns_gaps = path_output_directory + ClusterID + "_filter_dm_Ns_gaps.fas"
				retrieve_fasta_dict(list_record_ids,path_new_fasta,fasta_filtered_by_dist_matrix_Ns_gaps)
				#Third alignment
				translatorX_prefix_filtered_dist_Ns_gaps = ClusterID + "_filter_dm_Ns_gaps"
				output_trimal_dna,list_records_dna = wraper_translatorx_trimal(fasta_filtered_by_dist_matrix_Ns_gaps,translatorX_prefix_filtered_dist_Ns_gaps,path_output_directory)
				
				#Save the final filtered fasta file
				final_filtered_fasta = path_output_directory + ClusterID + "_filtered.fas"
				shutil.move(output_trimal_dna,final_filtered_fasta)
				#Remove intermediate files
				if cleanup == "yes":
					cleaning_intermediate_files(path_output_directory,ClusterID)
				logging.info(ClusterID + ": Has been successefully filtered.")
		else:
			logging.info(ClusterID + ": All sequences have been removed due to they are not multiple of 3 or unexpected stop codons")
			return

'''
Main
'''

#Avoid overwrite previous result and log files.
assert cleanup in ["yes","no"], "WARNING: The option cleanup must be yes or not."

assert not os.path.exists(log_file), "WARNING: The file " + log_file + " already exist. Please rename or move the previous log file to prevent information loss."
#assert not os.path.exists(out_file), "WARNING: The file " + out_file + " already exist. Please rename or move the previous results to prevent information loss."

#Create the log file.
logging.basicConfig(filename=log_file, filemode='w', format='%(levelname)s:%(message)s', level=logging.DEBUG)

#Run the program
cleaning(fasta_file,out_dir,distance_threshold,proportion_Ns_gaps,length_filter,cleanup)
