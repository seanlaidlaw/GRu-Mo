#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sqlite3


# Filepaths
fasta_path = "example_data/small/myfasta.fasta"
fastq_path = "example_data/small/myfastq_orig.fastq"
database_path = 'database.sqlite'

# Hyperparameters
kmer_size = 24
max_gap_size = 7
penalty_cutoff = 25
mismatch_penalty = 2
gap_opening_penalty = 1
gap_continuing_penalty = 1

os.system("rm 'database.sqlite'")
sqlite_connection = sqlite3.connect(database_path)


def calculate_penalty(sequenceA, sequenceB):
	current_score = 0
	lenA = len(sequenceA)
	lenB = len(sequenceB)
	smallest_len = lenA
	if lenB < lenA:
		smallest_len = lenB

	for i in range(1, smallest_len):
		a = sequenceA[i]
		b = sequenceB[i]
		if a == "-":
			if sequenceA[i-1] == "-":
				current_score += gap_continuing_penalty
			else:
				current_score += gap_opening_penalty

		elif b == "-":
			if sequenceB[i-1] == "-":
				current_score += gap_continuing_penalty
			else:
				current_score += gap_opening_penalty

		elif a != b:
			current_score += mismatch_penalty
	return current_score


def is_offset(fastxType):

	if fastxType == "fasta":
		fastxLine = fastaLine
		inverse_fastxLine = fastqLine
		fastx_rel_pos = fasta_rel_pos
		inverse_fastx_rel_pos = fastq_rel_pos

	elif fastxType == "fastq":
		fastxLine = fastqLine
		inverse_fastxLine = fastaLine
		fastx_rel_pos = fastq_rel_pos
		inverse_fastx_rel_pos = fasta_rel_pos

	else:
		print("Incorrect Argument passed to is_offset()")
		return False

	for gapSize in range(1, max_gap_size + 1):
		offset_gap = offset_gap_exists(gapSize, fastxLine, fastx_rel_pos, inverse_fastxLine, inverse_fastx_rel_pos)
		if offset_gap:  # is false if no gaps gives a better score
			return gapSize
	else:
		return False


def offset_gap_exists(gapSize, sequenceA, seqA_rel_pos, sequenceB, seqB_rel_pos):
	try:
		sequenceB[seqB_rel_pos + gapSize + 5]
		sequenceA[seqA_rel_pos + 5]
	except IndexError:
		return False

	if sequenceB[seqB_rel_pos + gapSize] == sequenceA[seqA_rel_pos]:
		# verify its a gap and not just 1/4 chance by testing next 5 nucleotides for offset
		if sequenceB[seqB_rel_pos + gapSize + 1] == sequenceA[seqA_rel_pos + 1]:
			if sequenceB[seqB_rel_pos + gapSize + 2] == sequenceA[seqA_rel_pos + 2]:
				if sequenceB[seqB_rel_pos + gapSize + 3] == sequenceA[seqA_rel_pos + 3]:
					if sequenceB[seqB_rel_pos + gapSize + 4] == sequenceA[seqA_rel_pos + 4]:
						if sequenceB[seqB_rel_pos + gapSize + 5] == sequenceA[seqA_rel_pos + 5]:
							return True

					else:
						return False
				else:
					return False
			else:
				return False
		else:
			return False


# Create database schema
db_cursor = sqlite_connection.cursor()
db_cursor.execute("CREATE TABLE IF NOT EXISTS ref_table(\
			  id INTEGER PRIMARY KEY  AUTOINCREMENT,\
			  quant INTEGER,\
			  header TEXT,\
			  transcripts VARCHAR(80));")


# Parse fasta file into SQLite3 database
with open(fasta_path) as myfasta:
	current_header = ""
	for line in myfasta:
		if line.startswith(">"):
			current_header = line.strip(">").strip("\n")
		else:
			if line is not " " and line is not "":
				line = line.strip("\n")
				db_cursor.execute(
					"INSERT INTO ref_table (transcripts, quant, header) VALUES ('" + line + "',0,'" + current_header + "');")

# Commit changes to database
sqlite_connection.commit()

# Index databse to allow Log n lookups
db_cursor.execute("CREATE INDEX index_name ON ref_table (transcripts);")

# Align fastq to database, read by read
with open(fastq_path) as myfastq:
	line_counter = -2
	for fastqLine in myfastq:
		line_counter += 1

		if line_counter % 4 == 0:  # to only read lines with DNA sequences
			fastqLine = fastqLine.strip("\n")
			kmer = fastqLine[0:kmer_size]

			# Get id of database row where kmer exists (less expensive option because O(log n) )
			fastaLine_ids = db_cursor.execute(
				"SELECT id FROM ref_table WHERE transcripts LIKE '" + kmer + "%';").fetchall()
			if len(fastaLine_ids) == 0:
				# using more expensive O(n) kmer search
				fastaLine_ids = db_cursor.execute(
					"SELECT id FROM ref_table WHERE transcripts LIKE '%" + kmer + "%';").fetchall()

			for fastaLine_id in fastaLine_ids:
				fastaLine_id = fastaLine_id[0]  # extract tuple from list
				fastaLine_seq = \
					db_cursor.execute(
						"SELECT transcripts FROM ref_table WHERE id=" + str(fastaLine_id) + ";").fetchall()[0]
				fastaLine_seq = fastaLine_seq[0]
				fastaLine = fastaLine_seq
				align_start = fastaLine.find(kmer)
				kmerend = align_start + kmer_size
				matched_seq_fasta = kmer
				matched_seq_fastq = kmer
				read_length = len(fastqLine)
				fastaLine = fastaLine[align_start:]

				for i in range(kmer_size, read_length):
					fastq_rel_pos = i  # starting point is just after end of kmer, end is at length of read
					fasta_rel_pos = fastq_rel_pos + align_start  # fasta position is based on where kmer was found

					fastaLine_id_plus1 = db_cursor.execute(
						"SELECT transcripts FROM ref_table where id=" + str(fastaLine_id + 1) + ";").fetchall()

					try:
						fastaLine[fasta_rel_pos]
					except:
						# if line from fasta is too small to compare whole read then get the following line from database
						fastaLine_id_plus1 = db_cursor.execute(
							"SELECT transcripts FROM ref_table where id=" + str(fastaLine_id + 1) + ";").fetchall()
						if fastaLine_id_plus1:
							fastaLine_id_plus1 = fastaLine_id_plus1[0]
							fastaLine_id_plus1 = fastaLine_id_plus1[0]
							fastaLine = str(fastaLine) + str(fastaLine_id_plus1)

					if fastqLine[fastq_rel_pos] == fastaLine[fasta_rel_pos]:
						matched_seq_fastq += fastaLine[fastq_rel_pos]
						matched_seq_fasta += fastaLine[fastq_rel_pos]

					elif is_offset("fasta"):
						gapSize = is_offset("fasta")  # returns size of gap
						fastaLine = fastaLine[0:fasta_rel_pos] + ("-" * gapSize) + fastaLine[fasta_rel_pos:]
						matched_seq_fasta += "-"
						matched_seq_fastq += fastqLine[fastq_rel_pos]


					elif is_offset("fastq"):
						gapSize = is_offset("fastq")  # returns size of gap
						fastqLine = fastqLine[0:fastq_rel_pos] + ("-" * gapSize) + fastqLine[fastq_rel_pos:]
						matched_seq_fastq += "-"
						matched_seq_fasta += fastaLine[fasta_rel_pos]

					else:  # mismatch
						matched_seq_fasta += fastaLine[fasta_rel_pos]
						matched_seq_fastq += fastqLine[fastq_rel_pos]

				my_penalty = calculate_penalty(matched_seq_fasta, matched_seq_fastq)
				diff_count = 0

				if my_penalty < penalty_cutoff:
					print("\nPenalty : " + str(my_penalty))
					print(matched_seq_fastq)
					print(matched_seq_fasta)
					db_cursor.execute("UPDATE ref_table SET quant = quant + 1 WHERE id=" + str(fastaLine_id) + ";")

sqlite_connection.commit()

quant_scores = db_cursor.execute("SELECT header,SUM(quant) FROM ref_table GROUP BY header;").fetchall()
print("\n Scores:")
for score in quant_scores:
	print(score)

sqlite_connection.close()
