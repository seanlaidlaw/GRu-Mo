#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sqlite3 # handle sqlite db
import sqlalchemy # handle sqlite db

os.system("rm 'database.sqlite'")
sqlite_connection = sqlite3.connect('database.sqlite')

kmersize = 14
maxgapsize = 4


def fastq_is_offset():
	for gapSize in range(1, maxgapsize + 1):
		if fastaLine[fasta_rel_pos + gapSize] == fastqLine[fastq_rel_pos]:
			# verify its a gap and not just 1/4 chance by testing next 4 nucleotides for decalage
			if fastaLine[fasta_rel_pos + gapSize + 1] == fastqLine[fastq_rel_pos + 1]:
				if fastaLine[fasta_rel_pos + gapSize + 2] == fastqLine[fastq_rel_pos + 2]:
					if fastaLine[fasta_rel_pos + gapSize + 3] == fastqLine[fastq_rel_pos + 3]:
						if fastaLine[fasta_rel_pos + gapSize + 4] == fastqLine[fastq_rel_pos + 4]:
							return gapSize

# Same function as above but checks fastaLine for offset
def fasta_is_offset():
	for gapSize in range(1, maxgapsize + 1):
		if fastqLine[fastq_rel_pos + gapSize] == fastaLine[fasta_rel_pos]:
			# verify its a gap and not just 1/4 chance by testing next 4 nucleotides for decalage
			if fastqLine[fastq_rel_pos + gapSize + 1] == fastaLine[fasta_rel_pos + 1]:
				if fastqLine[fastq_rel_pos + gapSize + 2] == fastaLine[fasta_rel_pos + 2]:
					if fastqLine[fastq_rel_pos + gapSize + 3] == fastaLine[fasta_rel_pos + 3]:
						if fastqLine[fastq_rel_pos + gapSize + 4] == fastaLine[fasta_rel_pos + 4]:
							return gapSize

# fastaLine = "GCAATGGGGCCCAACCCTTGGAGGCACTGCCCTTGGAGGC" + "TAT" + "AACGACCCGAAAATCTAGAACAGAAACCTAA"
# fastqLine =               "CCCTTGGAGGCACTGCCCTTGGAGGC" + ""    + "AACGACCCCGAAAATCAGTAACA"

db_cursor = sqlite_connection.cursor()
db_cursor.execute("CREATE TABLE ref_table(\
				  ID INTEGER PRIMARY KEY  AUTOINCREMENT,\
				  transcripts VARCHAR(80));")


with open("myfasta.fasta") as myfasta:
	for line in myfasta:
		if not line.startswith(">"):
			if line is not " " and line is not "":
				db_cursor.execute("INSERT INTO ref_table (transcripts) VALUES ('" + line.strip("\n") + "');")

db_cursor.execute("CREATE INDEX index_name ON ref_table (transcripts);")
sqlite_connection.commit()

with open("myfastq.fastq") as myfastq:
	for fastqLine in myfastq:
		kmer = fastqLine[0:kmersize]
		fastaLine_id = db_cursor.execute("SELECT id FROM ref_table where transcripts like '%" + kmer + "%';").fetchone()
		fastaLine_id = fastaLine_id[0]
		fastaLinelist = db_cursor.execute("SELECT transcripts FROM ref_table where id=" + str(fastaLine_id) + ";").fetchall()[0]
		fastaLinelist = fastaLinelist[0]
		fastaLinelist2 = db_cursor.execute("SELECT transcripts FROM ref_table where id=" + str(fastaLine_id + 1) + ";").fetchall()
		if fastaLinelist2:
			fastaLinelist2 = fastaLinelist2[0]
			fastaLinelist2 = fastaLinelist2[0]
			fastaLinelist = fastaLinelist + fastaLinelist2
		if len(fastaLinelist) == 1:
			fastaLine = fastaLinelist
			align_start = fastaLine.find(kmer)
			kmerend = align_start + kmersize
			matched_seq_fasta = kmer
			matched_seq_fastq = kmer
			print(matched_seq_fastq )
			read_length = len(fastqLine)

			# for i in range(7, 37):
			for i in range(kmersize, read_length):
				fastq_rel_pos = i  # starting point is just after end of kmer, end is at length of read
				fasta_rel_pos = fastq_rel_pos + align_start  # fasta position is based on where kmer was found
				print(fastaLine)
				print(fasta_rel_pos)
				print(fastqLine)
				print(fastq_rel_pos)

				if fastqLine[fastq_rel_pos] == fastaLine[fasta_rel_pos]:
					matched_seq_fastq += fastaLine[fastq_rel_pos]
					matched_seq_fasta += fastaLine[fastq_rel_pos]

				elif fasta_is_offset():
					gapSize = fasta_is_offset()  # returns size of gap
					fastaLine = fastaLine[0:fasta_rel_pos] + ("-" * gapSize) + fastaLine[fasta_rel_pos:]
					matched_seq_fasta += "-"
					matched_seq_fastq += fastqLine[fastq_rel_pos]

				elif fastq_is_offset():
					gapSize = fastq_is_offset()  # returns size of gap
					fastqLine = fastqLine[0:fastq_rel_pos] + ("-" * gapSize) + fastqLine[fastq_rel_pos:]
					matched_seq_fastq += "-"
					matched_seq_fasta += fastaLine[fasta_rel_pos]

				else:  # mismatch
					matched_seq_fasta += fastaLine[fasta_rel_pos]
					matched_seq_fastq += fastqLine[fastq_rel_pos]

			print(matched_seq_fastq)
			print(matched_seq_fasta)
sqlite_connection.close()
