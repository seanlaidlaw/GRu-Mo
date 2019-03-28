#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sqlite3

os.system("rm 'database.sqlite'")
sqlite_connection = sqlite3.connect('database.sqlite')

fasta_path = "example_data/small/myfasta.fasta"
fastq_path = "example_data/small/myfastq.fastq"

kmersize = 24
maxgapsize = 5
similarity_score_cutoff = 80


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

	for gapSize in range(1, maxgapsize + 1):
		try:
			if inverse_fastxLine[inverse_fastx_rel_pos + gapSize] == fastxLine[fastx_rel_pos]:
				# verify its a gap and not just 1/4 chance by testing next 4 nucleotides for decalage
				if inverse_fastxLine[inverse_fastx_rel_pos + gapSize + 1] == fastxLine[fastx_rel_pos + 1]:
					if inverse_fastxLine[inverse_fastx_rel_pos + gapSize + 2] == fastxLine[fastx_rel_pos + 2]:
						if inverse_fastxLine[inverse_fastx_rel_pos + gapSize + 3] == fastxLine[fastx_rel_pos + 3]:
							if inverse_fastxLine[inverse_fastx_rel_pos + gapSize + 4] == fastxLine[fastx_rel_pos + 4]:
								return gapSize
		except IndexError:
			return False


db_cursor = sqlite_connection.cursor()
db_cursor.execute("CREATE TABLE ref_table(\
			  ID INTEGER PRIMARY KEY  AUTOINCREMENT,\
			  quant INTEGER,\
			  header TEXT,\
			  transcripts VARCHAR(80));")

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

db_cursor.execute("CREATE INDEX index_name ON ref_table (transcripts);")
sqlite_connection.commit()
print("> parsed .fasta added to DB")

with open(fastq_path) as myfastq:
	line_counter = -2
	for fastqLine in myfastq:
		line_counter += 1
		if line_counter % 4 == 0:  # to only read lines with DNA sequences
			fastqLine = fastqLine.strip("\n")
			kmer = fastqLine[0:kmersize]
			# Get id of database row where kmer exists (less expensive option because O(log n) )
			fastaLine_ids = db_cursor.execute(
				"SELECT id FROM ref_table where transcripts like '" + kmer + "%';").fetchall()
			if len(fastaLine_ids) == 0:
				# using more expensive O(N) kmer search
				fastaLine_ids = db_cursor.execute(
					"SELECT id FROM ref_table where transcripts like '%" + kmer + "%';").fetchall()

			for fastaLine_id in fastaLine_ids:
				fastaLine_id = fastaLine_id[0]  # extract tuple from list
				fastaLine_seq = \
					db_cursor.execute(
						"SELECT transcripts FROM ref_table where id=" + str(fastaLine_id) + ";").fetchall()[0]
				fastaLine_seq = fastaLine_seq[0]
				fastaLine = fastaLine_seq
				align_start = fastaLine.find(kmer)
				kmerend = align_start + kmersize
				matched_seq_fasta = kmer
				matched_seq_fastq = kmer
				read_length = len(fastqLine)

				# for i in range(7, 37):
				for i in range(kmersize, read_length):
					fastq_rel_pos = i  # starting point is just after end of kmer, end is at length of read
					fasta_rel_pos = fastq_rel_pos + align_start  # fasta position is based on where kmer was found
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

				diff_count = 0
				for a, b in zip(matched_seq_fastq, matched_seq_fasta):
					if a != b or a == "-" or b == "-":
						diff_count += 1
				similarity_score = 100 * (read_length - diff_count) / read_length

				print("\nSimilarity : " + str(similarity_score) + "%")
				print(matched_seq_fastq)
				print(matched_seq_fasta)
				if similarity_score > similarity_score_cutoff:
					db_cursor.execute("UPDATE ref_table SET quant = quant + 1 WHERE id=" + str(fastaLine_id) + ";")

sqlite_connection.commit()

quant_scores = db_cursor.execute("SELECT header,SUM(quant) FROM ref_table GROUP BY header;").fetchall()
for score in quant_scores:
	print(score)

sqlite_connection.close()
