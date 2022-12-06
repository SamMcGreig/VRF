import sys
sys.path.insert(0, '/home/smcgreig/Scripts/Angua_test/')

import os
import subprocess
import Bio
import pysam
import argparse
import statistics
import pandas as pd
import re
import math
from collections import OrderedDict
from Angua import looper, create_dirs, document_env
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq

def main():

	# Find Variable Region script version
	fvr_version = 1

	# Input arguments
	options = parse_arguments()

	if(sys.argv[1] == "main"):

		# Directory list to create
		dirs = ["1.MetaEuk",
				"2.Bed_out",
				"3.Genes_out",
				"4.Proteinortho",
				"5.Proteinortho_groups",
				"6.Pseudo_reference",
				"7.GATK",
				"8.Variable_regions"
			]

		if options.create_dirs == "Y":
			create_dirs(options.output, dirs)

		if options.metaeuk == "Y":
			run_metaeuk(options.input_assemblies, options.input_ref_proteins, f"{options.output}/1.MetaEuk/", options.metaeuk_threads)	

		if options.bedtools == "Y":
			get_gene_coords(f"{options.output}/1.MetaEuk/", "results.fas", f"{options.output}/2.Bed_out", options.bedtools_delim)

			extract_genes(f"{options.output}/2.Bed_out/", options.input_assemblies, "_gapClosed.fa", f"{options.output}/3.Genes_out/")

		if options.proteinortho == "Y":
			run_proteinortho(f"{options.output}/3.Genes_out/", 
				options.proteinortho_threads, 
				options.name, 
				f"{options.output}/4.Proteinortho/", 
				options.genome_count, 
				options.proteinortho_similarity_cutoff, 
				f"{options.output}/3.Genes_out/", 
				f"{options.output}/5.Proteinortho_groups/")

			create_pseudo_ref(f"{options.output}/5.Proteinortho_groups/", f"{options.output}/6.Pseudo_reference/")

		if options.gatk == "Y":

			# Create directories
			dirs = ["1.bwa",
					"2.variant_calls",
					"3.extract_SNP_INDEL",
					"4.filter_SNP_INDEL",
					"5.select_variants",
					"6.recalibration",
					"7.variant_calls",
					"8.extract_SNP_INDEL",
					"9.filter_SNP_INDEL",
					"10.phased_SNP"
					]

			# Define SNP and INDEL dicts
			snp_dict = {"QD": f" < {options.snp_QD}",
						"FS": f" > {options.snp_FS}",
						"MQ": f"< {options.snp_MQ}",
						"SOR": f" > {options.snp_SOR}",
						"MQRankSum": f" < {options.snp_MQRankSum}",
						"ReadPosRankSum": f" < {options.snp_ReadPosRankSum}" 
						}
			indel_dict = {	"QD": f" < {options.indel_QD}",
							"FS": f" > {options.indel_FS}",
							"SOR": f" > {options.indel_SOR}"
						}

			# Preprocessing
			create_dirs(f"{options.output}/7.GATK/", dirs)
			run_bwa(options.input_raw, f"{options.output}/7.GATK/{dirs[0]}/", f"{options.output}/6.Pseudo_reference/pseudo_reference.fasta", options.bwa_threads)
			mark_duplicates(f"{options.output}/7.GATK/{dirs[0]}/")
			merge_bam(f"{options.output}/7.GATK/{dirs[0]}/")

			# GATK workflow
			gatk_run = GATK(f"{options.output}/6.Pseudo_reference/pseudo_reference.fasta", "samples.bam", f"{options.output}/7.GATK/", dirs)
			gatk_run.create_indexes()

			# Round 1 variant calling
			gatk_run.call_variants("raw")
			gatk_run.select_variants("SNP", "raw")
			gatk_run.select_variants("INDEL", "raw")
			gatk_run.filter_variants(snp_dict, "SNP", "raw")
			gatk_run.filter_variants(indel_dict, "INDEL", "raw")
			gatk_run.exclude_variants("SNP")
			gatk_run.exclude_variants("INDEL")

			# Recalibration
			gatk_run.recalibrate()

			# Round 2 variant calling
			gatk_run.call_variants("recalibrated")
			gatk_run.select_variants("SNP", "recalibrated")
			gatk_run.select_variants("INDEL", "recalibrated")
			gatk_run.filter_variants(snp_dict, "SNP", "recalibrated")
			gatk_run.filter_variants(indel_dict, "INDEL", "recalibrated")

			# Phase variants
			gatk_run.phase()

		if options.find_variable_regions == "Y":
			filter_snps(f"{options.output}/7.GATK/10.phased_SNP/phased_SNP.vcf", options.min_depth, options.marker_distance, options.variable_positions, f"{options.output}/8.Variable_regions/")
			select_informative_markers(f"{options.output}/8.Variable_regions/potential_markers.tsv", f"{options.output}/8.Variable_regions/marker_information.tsv")

		if options.doc_env == "Y":
			document_env("FindVariableRegions-main", fvr_version, options.output, options)

	elif(sys.argv[1] == "primer-design"):

		# Directory list to create
		dirs = ["Primer-design",
				"Primer-design/1.Bed_out",
				"Primer-design/2.Regions_out",
				"Primer-design/3.Primer3_input",
				"Primer-design/4.Primer3_output"
				]

		if options.create_dirs == "Y":
			create_dirs(options.output, dirs)

		if options.bedtools == "Y":
			create_bedfile(options.input_markers, f"{options.output}/Primer-design/1.Bed_out/markers.bed", f"{options.output}/Primer-design/1.Bed_out/unbuffered_markers.bed", 300)
			subprocess.call(f"bedtools getfasta -s -bed {options.output}/Primer-design/1.Bed_out/markers.bed -fi {options.input_reference} -fo {options.output}/Primer-design/2.Regions_out/markers.fasta", shell = True)

		if options.primer3 == "Y":
			generate_primer3_input(f"{options.output}/Primer-design/2.Regions_out/markers.fasta", f"{options.output}/Primer-design/3.Primer3_input/markers.txt", f"{options.output}/Primer-design/1.Bed_out/unbuffered_markers.bed")
			subprocess.call(f"primer3_core {options.output}/Primer-design/3.Primer3_input/markers.txt > {options.output}/Primer-design/4.Primer3_output/markers.txt", shell = True)
			collate_primers(f"{options.output}/Primer-design/4.Primer3_output/markers.txt", f"{options.output}/Primer-design/markers.tsv")

		if options.doc_env == "Y":
			document_env("FindVariableRegions-primer-design", fvr_version, options.output, options)

################################################################################
def run_metaeuk(assembly_dir, ref_proteins, output_dir, threads):
	input_files = looper(assembly_dir, "single", ".fa")

	for file in input_files:
		subprocess.call(f"metaeuk easy-predict {assembly_dir}/{file} {ref_proteins} {output_dir}/{file.split('.')[0]}.results {output_dir}/tmp --threads {threads}", shell = True)

	print("MetaEuk complete.")

################################################################################
def get_gene_coords(input_dir, qualifiier, output_dir, delim):
	input_files = looper(input_dir, "qualified", qualifiier)

	for file in input_files:
		with open(f"{input_dir}/{file}") as metaeuk_in:
			with open(f"{output_dir}/{file.split(delim)[0]}.bed", "w") as output_bedfile:
				for line in metaeuk_in:
					if ">" in line:
						line = line.strip()
						linesplit = line.split("|")
						protein = linesplit[0]
						seqID = linesplit[1]
						start = linesplit[6]
						end = linesplit[7]

						if abs(int(end) - int(start)) >= 300:
							output_bedfile.write(f"{seqID}	{start}	{end}\n")

	print("Gene coordinate extraction complete.")

################################################################################
def extract_genes(input_dir, assembly_dir, assembly_suffix, output_dir):
	input_files = looper(input_dir, "single", "bed")

	for file in input_files:
		sample_name = file.split(".")[0]
		subprocess.call(f"bedtools getfasta -s -bed {input_dir}/{file} -fi {assembly_dir}/{sample_name}{assembly_suffix} -fo {output_dir}/{sample_name}.fasta", shell = True)

	print("Gene extraction complete.")

################################################################################
def run_proteinortho(input_dir, threads, name, output_dir, genome_count, similarity_cutoff, gene_out_dir, output_dir2):
	subprocess.call(f"proteinortho {input_dir}/*.fasta -p=blastn -cpus={threads} -clean -project={name}", shell = True)
	subprocess.call(f"mv {name}* {output_dir}/", shell = True)

	with open(f"{output_dir}/{name}.proteinortho.tsv") as ortho_out:
		with open(f"{output_dir}/{name}.proteinortho.filtered.tsv", "w") as proteinortho_filtered_out:
			header = next(ortho_out)
			proteinortho_filtered_out.write(header)
			for line in ortho_out:
				linestrip = line.split("\t")

				if linestrip[0] == genome_count:
					if linestrip[1] == genome_count:
						if float(linestrip[2]) >= float(similarity_cutoff):
							proteinortho_filtered_out.write(line.replace(")", "))"))

	subprocess.call(f"proteinortho_grab_proteins.pl -tofiles {output_dir}/{name}.proteinortho.filtered.tsv -source {gene_out_dir}/*", shell = True)
	subprocess.call(f"for file in {name}.proteinortho.filtered.tsv.OrthoGroup*; do mv $file {output_dir2}/$file; done", shell = True)

	print("Proteinortho complete.")

################################################################################
def create_pseudo_ref(input_dir, output_dir):
	input_files = looper(input_dir, "qualified", ".fasta")

	with open(f"{output_dir}/pseudo_reference.fasta", "a") as pseudo_ref_out:
		for file in input_files:
			with open(f"{input_dir}/{file}") as in_fasta:
				for record in SeqIO.parse(in_fasta, "fasta"):
					pseudo_ref_out.write(f">{record.id[:-2]}\n")
					pseudo_ref_out.write(f"{record.seq}\n")
					break

################################################################################
def run_bwa(input_dir, output_dir, reference, threads):
	input_files = looper(input_dir, "qualified", "_R1_001.fastq.gz")

	subprocess.call(f"bwa index {reference}", shell = True)

	for file in input_files:
		r1 = f"{input_dir}/{file}"
		r2 = f"{input_dir}/{file.replace('_R1', '_R2')}"
		sample_name = file.split("_S")[0]
		outfile = f"{output_dir}/{sample_name}.aligned_reads.sam"

		subprocess.call(f"bwa mem -R \'@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPM:NOVASEQ\' -K 100000000 -Y -t {threads} {reference} {r1} {r2} > {outfile}", shell = True)

	print("BWA completed.")

################################################################################
def mark_duplicates(input_dir):
	input_files = looper(input_dir, "qualified", ".sam")

	for file in input_files:
		infile = f"{input_dir}/{file}"
		sample_name = file.split(".")[0]
		outfile = f"{input_dir}/{sample_name}.sorted_dedup_reads.bam"

		subprocess.call(f"gatk MarkDuplicatesSpark -I {infile} -O {outfile}", shell = True)

		# Clean up
		os.remove(f"{input_dir}/{file}")

	print("Duplicates marked and sorted.")

################################################################################
def merge_bam(input_dir):
	input_files = looper(input_dir, "qualified", ".bam")

	input_string = ""
	for file in input_files:
		if file.endswith(".bam"):
			infile = f"{input_dir}/{file}"

			input_string = f"{input_string}I={infile} "

	subprocess.call(f"picard MergeSamFiles {input_string} O={input_dir}/samples.bam", shell = True)
	subprocess.call(f"samtools index {input_dir}/samples.bam", shell = True)

	# Clean up
	for file in input_files:
		infile = f"{input_dir}/{file}"

		os.remove(f"{input_dir}/{file}")

	print("bam files merged.")

################################################################################
class GATK:
	def __init__(self, reference, bwa_filename, gatk_root_dir, dirs):
		self.reference = reference
		self.bwa_filename = bwa_filename
		self.dirs = [gatk_root_dir + directory for directory in dirs]

	############################################################################
	def create_indexes(self):
		subprocess.call(f"samtools faidx {self.reference}", shell = True)
		subprocess.call(f"gatk CreateSequenceDictionary REFERENCE={self.reference} OUTPUT={self.reference.split('.fasta')[0]}.dict", shell = True)

	############################################################################
	def call_variants(self, mode):
		if mode == "raw":
			subprocess.call(f"gatk HaplotypeCaller -R {self.reference} -I {self.dirs[0]}/{self.bwa_filename} -O {self.dirs[1]}/{mode}_variants.vcf", shell = True)
		elif mode == "recalibrated":
			subprocess.call(f"gatk HaplotypeCaller -R {self.reference} -I {self.dirs[5]}/recallibrated_reads.bam -O {self.dirs[6]}/{mode}_variants.vcf", shell = True)

	############################################################################
	def select_variants(self, variant_type, mode):
		if mode == "raw":
			subprocess.call(f"gatk SelectVariants -R {self.reference} -V {self.dirs[1]}/{mode}_variants.vcf -select-type {variant_type} -O {self.dirs[2]}/{mode}_{variant_type}.vcf", shell = True)
		elif mode == "recalibrated":
			subprocess.call(f"gatk SelectVariants -R {self.reference} -V {self.dirs[6]}/{mode}_variants.vcf -select-type {variant_type} -O {self.dirs[7]}/{mode}_{variant_type}.vcf", shell = True)

	############################################################################
	def filter_variants(self, filter_dict, variant_type, mode):
		filter_string = ""
		for filter_param in filter_dict:
			filter_string = f"{filter_string}-filter-name \"{filter_param}_filter\" -filter \"{filter_param}{filter_dict[filter_param]}\" "

		if mode == "raw":
			subprocess.call(f"gatk VariantFiltration -R {self.reference} -V {self.dirs[2]}/{mode}_{variant_type}.vcf -O {self.dirs[3]}/filtered_{variant_type}.vcf {filter_string}", shell = True)
		elif mode == "recalibrated":
			subprocess.call(f"gatk VariantFiltration -R {self.reference} -V {self.dirs[7]}/{mode}_{variant_type}.vcf -O {self.dirs[8]}/filtered_{variant_type}.vcf {filter_string}", shell = True)

	############################################################################
	def exclude_variants(self, variant_type):
		subprocess.call(f"gatk SelectVariants --exclude-filtered -V {self.dirs[3]}/filtered_{variant_type}.vcf -O {self.dirs[4]}/bqsr_{variant_type}.vcf", shell = True)

	############################################################################
	def recalibrate(self):
		subprocess.call(f"gatk BaseRecalibrator -R {self.reference} -I {self.dirs[0]}/{self.bwa_filename} --known-sites {self.dirs[4]}/bqsr_SNP.vcf --known-sites {self.dirs[4]}/bqsr_INDEL.vcf -O {self.dirs[5]}/recallibrated_data.table", shell = True)
		subprocess.call(f"gatk ApplyBQSR -R {self.reference} -I {self.dirs[0]}/{self.bwa_filename} -bqsr {self.dirs[5]}/recallibrated_data.table -O {self.dirs[5]}/recallibrated_reads.bam", shell = True)

	############################################################################
	def phase(self):
		subprocess.call(f"whatshap phase -o {self.dirs[9]}/phased_SNP.vcf --reference {self.reference} {self.dirs[8]}/filtered_SNP.vcf {self.dirs[5]}/recallibrated_reads.bam", shell = True)

################################################################################
def filter_snps(snps, min_depth, marker_distance, variable_positions, output_dir):
	header = ["CONTIG"]

	snp_dict = OrderedDict()
	with open(snps) as snps_in:
		for line in snps_in:
			line = line.strip()
			# Remove unneeded headers
			if line.startswith("#CHROM"):
				linesplit = line.split("\t")
				del linesplit[0]
				del linesplit[1]
				del linesplit[3]
				del linesplit[4:6]
				del linesplit[3]
				header.extend(linesplit)

			# Create dictionary with lists
			elif not line.startswith("#"):
				linesplit = line.split("\t")

				# Create nested dictionary structure
				snp_dict.setdefault(linesplit[0], OrderedDict())
				snp_dict[linesplit[0]].setdefault(linesplit[1], [])

				# Collect static variables and append to the list
				ref, alt, filt = linesplit[3], linesplit[4], linesplit[6]
				snp_dict[linesplit[0]][linesplit[1]].extend([ref, alt, filt])

				# Collect sample variables and append to the list
				for i in range(9, len(linesplit)):
					genotype = linesplit[i].split(":")[0]
					depth = linesplit[i].split(":")[2]
					snp_dict[linesplit[0]][linesplit[1]].extend([genotype, depth])

	# Remove any sequences which have failed calls
	# Remove any sequences which have depth below X
	failed_calls = []
	low_depth = []
	for seq in snp_dict:
		for position in snp_dict[seq]:
			if snp_dict[seq][position][2] != "PASS":
				failed_calls.append(seq)
				break

			depth = []
			for value in snp_dict[seq][position][3:]:
				if "/" not in value and "|" not in value:
					if "." not in value:
						depth.append(int(value))

			if float(min_depth) > statistics.median(depth):
				low_depth.append(seq)
				break

	# Combine fails list and remove kets from dictionary
	seqs_to_remove = list(set(failed_calls + low_depth))
	for seq_to_remove in seqs_to_remove:
		snp_dict.pop(seq_to_remove)

	# Print snp_dict
	# print(snp_dict)

	# Find variable regions
	variable_position_dict = OrderedDict()
	for seq in snp_dict:
		snp_positions = [0]
		counter = 0
		for position in snp_dict[seq]:
			# If the first position is a 0, replace it with first SNP position we have
			if snp_positions[0] == 0:
				snp_positions = [int(position)]

			# If the first position is not a 0, then check to see how close the current position is to the previous position
			# If it is within the specified marker distance, append this to the list and add 1 to the counter
			# If this is condtion is not met, the counter and list of positions is reset
			elif snp_positions[0] != 0:
				if snp_positions[-1] + int(marker_distance) >= int(position):
					snp_positions.append(int(position))
					counter += 1
				elif int(position) > snp_positions[-1] + int(marker_distance):
					if counter >= int(variable_positions):
						# print(f"{seq}	{snp_positions}")
						variable_position_dict.setdefault(seq, [])
						for pos in snp_positions:
							variable_position_dict[seq].append(str(pos))
					snp_positions = [0]
					counter = 0

			# When the end of the cycle is detected, if the counter threshold is met then the result is reported
			if int(position) == int(list(snp_dict[seq].keys())[-1]):
				if counter >= int(variable_positions):
					# print(f"{seq}	{snp_positions}")
					variable_position_dict.setdefault(seq, [])
					for pos in snp_positions:
						variable_position_dict[seq].append(str(pos))
					snp_positions = [0]
					counter = 0

	# Report potential markers
	with open(f"{output_dir}/potential_markers.tsv", "w") as outfile:
		outfile.write(f"{'	'.join(header)}\n")
		for seq in variable_position_dict:
			for pos in variable_position_dict[seq]:
				# Collect static variables
				ref = snp_dict[seq][pos][0]
				alt = snp_dict[seq][pos][1]

				output = f"{seq}	{pos}	{ref}	{alt}"
				# Collect sample variables and add to output string
				for i in range(3, len(snp_dict[seq][pos])):
					if "/" in snp_dict[seq][pos][i] or "|" in snp_dict[seq][pos][i]:
						genotype = snp_dict[seq][pos][i].split(":")[0]
						output = f"{output}	{genotype}"
				output = f"{output}\n"
				outfile.write(output)

	print("SNPs filtered.")

################################################################################
def select_informative_markers(input_markers, output_information):

	# Read in as a pandas dataframe
	marker_data = pd.read_csv(input_markers, sep = "\t")

	# Obtain the sequence IDs and remove duplicates
	seq_ids = list(set(list(marker_data["CONTIG"])))

	# Informative markers dict
	informative_markers_dict = OrderedDict()
	
	# For each sequence ID, create a subset for SNPs associated with that ID
	for seq in seq_ids:
		subset = marker_data[marker_data["CONTIG"] == seq]

		# Create an alleles dictionary to store the various alleles observered
		alleles = OrderedDict()
		count = 1
		genotypes = []
		for col in subset.columns:
			if col not in ["CONTIG", "POS", "REF", "ALT"]:

				# Create a new dataframe andsplit the columns by genotype
				df = pd.DataFrame()
				df[["HT1", "HT2"]] = subset[col].str.split(r'\D', expand = True)
				ht1 = "".join(list(df["HT1"]))
				ht2 = "".join(list(df["HT2"]))
				
				# For each haplotype, see if it has been observed yet
				# If not, add it to the alleles dictionary
				for ht in [ht1, ht2]:
					if ht not in alleles:
						alleles.setdefault(ht, count)
						count += 1

				# Create a list of lists, storng each genotype
				genotypes.append(sorted([alleles[ht1], alleles[ht2]]))


		informative_markers_dict[seq] = len(list(map(list,set(map(tuple, genotypes)))))

	with open(output_information, "w") as outfile:
		outfile.write("SEQUENCE	GENOTYPES\n")
		for seq in informative_markers_dict:
			outfile.write(f"{seq}	{informative_markers_dict[seq]}\n")

################################################################################
def create_bedfile(input_markers, output_bed, unbuffered_bed, amplicon_size):

	with open(output_bed, "w") as marker_bed:
		with open(unbuffered_bed, "w") as unbuffered_marker_bed:
			with open(input_markers) as markers_in:
				next(markers_in)
				for line in markers_in:
					line = line.strip()
					contig = line.split("\t")[0]
					start = int(line.split("\t")[1])
					end = int(line.split("\t")[2])
					marker_buffer = math.floor((int(amplicon_size) - (end - start)) / 2)
					marker_bed.write(f"{contig}	{start - marker_buffer}	{end + marker_buffer}\n")
					unbuffered_marker_bed.write(f"{contig}	{start}	{end}\n")

################################################################################
def generate_primer3_input(input_markers, output_primer3, unbuffered_markers):

	# Unbuffered dict
	unbuffered_dict = {}
	with open(unbuffered_markers) as unbuffered:
		for line in unbuffered:
			line = line.strip()
			unbuffered_dict.setdefault(line.split("\t")[0], [])
			unbuffered_dict[line.split("\t")[0]].append(line.split("\t")[1])
			unbuffered_dict[line.split("\t")[0]].append(line.split("\t")[2])

	with open(output_primer3, "w") as primer_input:
			for record in SeqIO.parse(input_markers, "fasta"):
				seq_name = record.name[:-2]
				primer_input.write(f"SEQUENCE_ID={seq_name}\n")
				primer_input.write(f"SEQUENCE_TEMPLATE={record.seq}\n")
				

				if ":".join(seq_name.split(":")[:-1]) in unbuffered_dict:
					buffered_start = int(seq_name.split(":")[-1].split("-")[0])
					vr_start = int(unbuffered_dict[":".join(seq_name.split(":")[:-1])][0])
					vr_end = int(unbuffered_dict[":".join(seq_name.split(":")[:-1])][1])
					buffer_size = vr_start - buffered_start
					vr_size = vr_end - vr_start

					primer_input.write(f"SEQUENCE_EXCLUDED_REGION={buffer_size - 1},{vr_size + 1}\n")

				primer_input.write("=\n")
					
################################################################################
def collate_primers(input_markers, output_markers):

	with open(output_markers, "w") as primer_output:
		primer_output.write("SEQUENCE_ID	SEQUENCE_TEMPLATE	PRIMER_PAIR_PENALTY	PRIMER_LEFT_SEQUENCE	PRIMER_RIGHT_SEQUENCE(Reverse Complement)	PRIMER_LEFT	PRIMER_RIGHT	PRIMER_LEFT_TM	PRIMER_RIGHT_TM	PRIMER_LEFT_GC_PERCENT	PRIMER_RIGHT_GC_PERCENT\n")
		primer_details = []
		with open(input_markers) as primer_input:
			for line in primer_input:
				line = line.strip()
				stat = line.split("=")[0]
				value = line.split("=")[1]

				targets = [
							"SEQUENCE_ID",
							"SEQUENCE_TEMPLATE",
							"PRIMER_PAIR_0_PENALTY",
							"PRIMER_LEFT_0_SEQUENCE",
							"PRIMER_RIGHT_0_SEQUENCE",
							"PRIMER_LEFT_0",
							"PRIMER_RIGHT_0",
							"PRIMER_LEFT_0_TM",
							"PRIMER_RIGHT_0_TM",
							"PRIMER_LEFT_0_GC_PERCENT",
							"PRIMER_RIGHT_0_GC_PERCENT"
						]
				if(stat in targets):
					primer_details.append(value)

				if line.startswith("="):
					if len(primer_details) > 2:
						primer_output.write(f"{'	'.join(primer_details)}\n")
					else:
						primer_output.write(f"{'	'.join(primer_details)}	FAIL\n")
					primer_details = []
					
################################################################################
def parse_arguments():
	parser = argparse.ArgumentParser(prog = "Variable Region Finder", description = "Runs the Variable Region Finder pipeline.")
	subparsers = parser.add_subparsers(help = "sub-command help")
	main = subparsers.add_parser("main", help = "Runs the Variable Region Finder pipeline.")
	primer_design = subparsers.add_parser("primer-design", help = "Runs the primer-design pipeline.")

	### Main

	# Key arguments
	main.add_argument("--input_assemblies", help = "This is the location of the assembly data directory.", required = True)
	main.add_argument("--input_raw", help = "This is the location of the raw data directory. This should include a subdirectory for each sample's paired end data.", required = True)
	main.add_argument("--input_ref_proteins", help = "This is the location of the reference protein file.", required = True)
	main.add_argument("--output", help = "This is where the output data will be generated.", required = True)
	main.add_argument("--name", help = "This is the project name. Shoud be one word only (no spaces, underscores etc).", required = True, default = "MyProject")

	# Extra arguments, useful for if a specific job has failed and you don't want to start from scratch
	main.add_argument("--create_dirs", help = "Creates the directory structure. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--metaeuk", help = "Runs metaeuk. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--bedtools", help = "Runs bedtools to extract gene locations. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--proteinortho", help = "Runs proteinortho. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--gatk", help = "Runs gatk. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--find_variable_regions", help = "Finds variable regions from a vcf file. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--doc_env", help = "Creates environment and parameter details. Default Y.", choices = ["Y", "N"], default = "Y")

	# Tool specific parameters

	# MetaEuk
	main.add_argument("--metaeuk_threads", help = "The number of threads to run MetaEUK with. Default round(0.9*total threads).", default = round(os.cpu_count() * 0.9))

	# Bedtools
	main.add_argument("--bedtools_delim", help = "The separator between sample name and other elements in the filename.", required = True)

	# Proteinortho
	main.add_argument("--proteinortho_threads", help = "The number of threads to run ProteinOrtho with. Default round(0.9*total threads).", default = round(os.cpu_count() * 0.9))
	main.add_argument("--genome_count", help = "The number of assemblies/genomes a gene must be present in to be included.", required = True)
	main.add_argument("--proteinortho_similarity_cutoff", help = "The similarity cutoff for genes to be considered orthologues. Default 0.85", default = 0.85)

	# GATK
	main.add_argument("--bwa_threads", help = "The number of threads to run bwa with. Default round(0.9*total threads).", default = round(os.cpu_count() * 0.9))

	main.add_argument("--snp_QD", help = "QD filter. QD < X. Default 2.0", default = 2.0)
	main.add_argument("--snp_FS", help = "FS filter. FS > X. Default 60.0", default = 60.0)
	main.add_argument("--snp_MQ", help = "MQ filter. MQ < X. Default 40.0", default = 40.0)
	main.add_argument("--snp_SOR", help = "SOR filter. SOR > X. Default 3.0", default = 3.0)
	main.add_argument("--snp_MQRankSum", help = "MQRankSum filter. MQRankSum < -X. Default X", default = -12.5)
	main.add_argument("--snp_ReadPosRankSum", help = "ReadPosRankSum filter. ReadPosRankSum < X. Default -8.0", default = -8.0)

	main.add_argument("--indel_QD", help = "QD filter. QD < X. Default 2.0", default = 2.0)
	main.add_argument("--indel_FS", help = "FS filter. FS > X. Default 200.0", default = 200.0)
	main.add_argument("--indel_SOR", help = "SOR filter. SOR > X. Default 10.0", default = 10.0)

	# Find Variable Regions
	main.add_argument("--min_depth", help = "The minimum median depth of a variant. Default 40.0", default = 40.0)
	main.add_argument("--marker_distance", help = "The maximum distance allowed between 2 markers. Default 50", default = 50)
	main.add_argument("--variable_positions", help = "The minimum number of variable regiosn that need to be present in a variable region. Default 3", default = 3)

	### Primer Design

	# Key arguments
	primer_design.add_argument("--input_reference", help = "This is the location of the reference sequences.", required = True)
	primer_design.add_argument("--input_markers", help = "This is the location of the input marker file.", required = True)
	primer_design.add_argument("--output", help = "This is where the output data will be generated.", required = True)

	# Extra arguments, useful for if a specific job has failed and you don't want to start from scratch
	primer_design.add_argument("--create_dirs", help = "Creates the directory structure. Default Y.", choices = ["Y", "N"], default = "Y")
	primer_design.add_argument("--bedtools", help = "Runs bedtools to extract gene locations. Default Y.", choices = ["Y", "N"], default = "Y")
	primer_design.add_argument("--primer3", help = "Runs primer3 to create primer sets for the markers. Default Y.", choices = ["Y", "N"], default = "Y")

	primer_design.add_argument("--doc_env", help = "Creates environment and parameter details. Default Y.", choices = ["Y", "N"], default = "Y")


	return parser.parse_args()

################################################################################

if __name__ == '__main__':
	main()