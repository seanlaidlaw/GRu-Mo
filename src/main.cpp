#include <chrono>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <regex>
#include <sqlite3.h>
#include "cxxopts.hpp"
#include "database_management.h"

auto start = chrono::steady_clock::now();

const char *database_path = "transcript_reference_index.sqlite";
vector<double> similarity_scores;

bool skip_database_creation;
bool skip_quant;
int verbosity_level;
double kmer_size;
double max_gap_size;
double penalty_cutoff;
double mismatch_penalty;
double gap_opening_penalty;
double gap_continuing_penalty;
string fasta_path;
string fastq_path;

// this replaces the string * x function in python
string insert_char_x_times(char input_char, int times) {
	string new_string;
	for (int i = 1; i <= times; ++i) {
		new_string += input_char;
	}
	return new_string;
}

string reverseComplement(string forwardSequence) {
			replace(forwardSequence.begin(), forwardSequence.end(), 'A', 'Z'); // replace all 'x' to 'y'
			replace(forwardSequence.begin(), forwardSequence.end(), 'T', 'A'); // replace all 'x' to 'y'
			replace(forwardSequence.begin(), forwardSequence.end(), 'Z', 'T'); // replace all 'x' to 'y'
			replace(forwardSequence.begin(), forwardSequence.end(), 'C', 'Z'); // replace all 'x' to 'y'
			replace(forwardSequence.begin(), forwardSequence.end(), 'G', 'C'); // replace all 'x' to 'y'
			replace(forwardSequence.begin(), forwardSequence.end(), 'Z', 'G'); // replace all 'x' to 'y'
		return forwardSequence;
}

// Assign a penalty score to a given alignment based on program arguments
double calculate_penalty(string sequenceA, string sequenceB) {
	double current_score = 0;
	unsigned int lenA = sequenceA.size();
	unsigned int lenB = sequenceB.size();
	int smallest_len = lenA;
	if (lenB < lenA)
		smallest_len = lenB;

	for (int i = 0; i < smallest_len; ++i) {
		char a = sequenceA[i];
		char b = sequenceB[i];
		if (a == '-')
			if (sequenceA[i - 1] == '-')
				current_score += gap_continuing_penalty;
			else
				current_score += gap_opening_penalty;

		else if (b == '-')
			if (sequenceB[i - 1] == '-')
				current_score += gap_continuing_penalty;
			else
				current_score += gap_opening_penalty;

		else if (a != b)
			current_score += mismatch_penalty;
	}

	return current_score;
}


// Calculate objective similarity (not based on penalty variables).
// For every identical nucleotide 1 point is added to the similarity score
double calculate_similarity(string sequenceA, string sequenceB) {
	double matched_nucleotides = 0;
	unsigned int lenA = sequenceA.size();
	unsigned int lenB = sequenceB.size();
	int smallest_len = lenA;
	int longest_len = lenB;
	if (lenB < lenA) {
		smallest_len = lenB;
		longest_len = lenA;
	}

	for (int i = 0; i < smallest_len; ++i) {
		char a = sequenceA[i];
		char b = sequenceB[i];

		if (a != '-' and b != '-') {
			if (a == b) {
				matched_nucleotides += 1;
			}
		}
	}
	double similarity = 100 * (matched_nucleotides / longest_len);

	return similarity;
}


// if the next 5 nt of fastaLine are the same as the next+1 5 nt of fastqLine then return true as is offset
bool offset_gap_exists(int gapSize, string sequenceA, int seqA_rel_pos, string sequenceB, int seqB_rel_pos) {
	if (((seqB_rel_pos + gapSize + 5) > sequenceB.size()) or ((seqA_rel_pos + 5) > sequenceA.size())) {
		return false;
	}

	// verify its a gap and not just 1/4 chance by testing next 5 nucleotides for offset
	if (sequenceB[seqB_rel_pos + gapSize] != sequenceA[seqA_rel_pos]) {
		return false;
	}
	if (sequenceB[seqB_rel_pos + gapSize + 1] != sequenceA[seqA_rel_pos + 1]) {
		return false;
	}
	if (sequenceB[seqB_rel_pos + gapSize + 2] != sequenceA[seqA_rel_pos + 2]) {
		return false;
	}
	if (sequenceB[seqB_rel_pos + gapSize + 3] != sequenceA[seqA_rel_pos + 3]) {
		return false;
	}
	if (sequenceB[seqB_rel_pos + gapSize + 4] != sequenceA[seqA_rel_pos + 4]) {
		return false;
	}
	if (sequenceB[seqB_rel_pos + gapSize + 5] != sequenceA[seqA_rel_pos + 5]) {
		return false;
	}
	return true;
}

int is_offset(string sequenceA, int seqA_rel_pos, string sequenceB, int seqB_rel_pos) {
	for (int gapSize = 1; gapSize < max_gap_size; ++gapSize) {
		bool offset_gap = offset_gap_exists(gapSize, sequenceA, seqA_rel_pos, sequenceB, seqB_rel_pos);
		if (offset_gap) {
			// is false if no gaps gives a better score
			return gapSize;
		}
	}
	return 0;
}


int main(int argc, char *argv[]) {
	cxxopts::Options options("GRu-Mo", "Grouped Read unifier, mapping optimiser");
	options.add_options()
			("skip-db", "Skip fasta parsing and database construction", cxxopts::value<bool>()->default_value("false"))
			("skip-quant", "Skip quant, just build index database", cxxopts::value<bool>()->default_value("false"))
			("v,verbose",
			 "Set verbosity level, 0 is no output, 1 is CSV of results, 2 is quant scores, 3 is read alignment",
			 cxxopts::value<int>()->default_value("3"))
			("a,fasta", "reference transcriptome input", cxxopts::value<string>())
			("q,fastq", "reads input", cxxopts::value<string>())
			("k,kmer-size", "kmer size to search fasta", cxxopts::value<int>()->default_value("24"))
			("c,cutoff", "cutoff above which read with this penalty won't be quantified",
			 cxxopts::value<int>()->default_value("25"))
			("g,max-gap", "maximum gapsize to look for gaps", cxxopts::value<int>()->default_value("7"))
			("m,mismatch", "penalty for nucleotide mismatch", cxxopts::value<int>()->default_value("2"))
			("o,gap-open", "penalty for opening a gap", cxxopts::value<int>()->default_value("1"))
			("e,gap-extend", "penalty for extending a gap", cxxopts::value<int>()->default_value("1"));

	auto result = options.parse(argc, argv);

	if (result.count("help") || result.arguments().size() < 1) {
		std::cout << options.help() << std::endl;
		return 0;
	}

	if (!result.count("fasta")) {
		cout << "No fasta file given, please provide one" << endl;
		return 1;
	}

	if (!result.count("fastq")) {
		cout << "No fastq file given, please provide one" << endl;
		return 1;
	}

	::skip_quant = result["skip-quant"].as<bool>();
	::skip_database_creation = result["skip-db"].as<bool>();
	::verbosity_level = result["verbose"].as<int>();
	::fasta_path = result["fasta"].as<string>();
	::fastq_path = result["fastq"].as<string>();
	::kmer_size = result["kmer-size"].as<int>();
	::max_gap_size = result["max-gap"].as<int>();
	::penalty_cutoff = result["cutoff"].as<int>();
	::mismatch_penalty = result["mismatch"].as<int>();
	::gap_opening_penalty = result["gap-open"].as<int>();
	::gap_continuing_penalty = result["gap-extend"].as<int>();


	sqlite3 *db;
	string sql_command;

	// Set db to be the connection to the database
	sqlite3_open(database_path, &db);

	ifstream fasta_file(fasta_path);
	ifstream fastq_file(fastq_path);

	// Return error if provided files doesn't exist
	if (!fasta_file.good()) {
		cout << "Error: Provided fasta file does not exist";
		return 1;
	}

	if (!fastq_file.good()) {
		cout << "Error: Provided fastq file does not exist";
		return 1;
	}

	sql_command = "UPDATE ref_table SET quant = 0 WHERE quant!=0;";
	sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);

	if (::skip_database_creation == false) {
		// Delete previously made tables
		sql_command = "DROP TABLE ref_table;";
		sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);

		// Create database schema
		sql_command = "CREATE TABLE IF NOT EXISTS ref_table(id INTEGER PRIMARY KEY  AUTOINCREMENT, quant INTEGER, header TEXT, transcripts VARCHAR(80), hash INTEGER);";
		sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);

		sqlite3_exec(db, "PRAGMA cache_size=10000", NULL, NULL, NULL);
		sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);


		string line;
		string current_header;
		vector<unsigned long long> hashTable;
		while (getline(fasta_file, line)) {
			if (line[0] == '>') {
				line = line.substr(1, line.size());
				line = regex_replace(line, regex("^ "), "");
				line = regex_replace(line, regex(" +"), " ");
				current_header = line;

			} else {

				// TODO: change the i+=1 to higher nombre to increse effeciancy
				for (int i = 0; i < (line.size() - kmer_size); i += 3) {
					string kmer_to_hash = line.substr(i, kmer_size); // get each kmer from line
					size_t kmer_hash = hash<std::string>()(kmer_to_hash); // get hash for each kmer

//					cout << kmer_to_hash << endl;
//					cout << kmer_hash << endl;

					sql_command =
							"INSERT INTO ref_table (transcripts, quant, header, hash) VALUES ('" + line + "',0,'" +
							current_header + "', '" + to_string(kmer_hash) + "');";
					sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);

				}
			}
		}
		sqlite3_exec(db, "END TRANSACTION", NULL, NULL, NULL);

		// Index database for log n lookup times
		sql_command = "CREATE INDEX hash_index ON ref_table (hash);CREATE INDEX id_index ON ref_table (id);";
		sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);
	}

	if (::skip_quant == false) {
		int line_counter = 0;
		int read_counter = 0;
		int aligned_counter = 0;
		string fastqLine;
		while (getline(fastq_file, fastqLine)) {
			line_counter += 1;
			if (line_counter % 4 == 2) {

				// make vector of fastqLine and its reverse complement
				vector<string> fastqLines;
				fastqLines.push_back(fastqLine);
//				string revfastqLine = reverseComplement(fastqLine);
//				fastqLines.push_back(revfastqLine);

				for (string fastqLine : fastqLines) {
					read_counter += 1;
					for (int i = 0; i < (fastqLine.size() - kmer_size); i += 3) {
						string kmer = fastqLine.substr(i, kmer_size);
						unsigned long long kmer_hash = hash<string>()(kmer); // get hash for each possible fastq kmer

						// Get id of database row where kmer exists (less expensive option because O(log n) )
						sql_command = "SELECT id FROM ref_table WHERE hash=" + to_string(kmer_hash) + ";";
						Records fastaLine_ids = select_stmt(sql_command.c_str(), db);

						for (auto &fastaLine_id : fastaLine_ids) {
//								sql_command = "UPDATE ref_table SET quant = quant + 1 WHERE id=" + fastaLine_id[0] + ";";
//								sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);
							sql_command = "SELECT transcripts FROM ref_table WHERE id=" + fastaLine_id[0] + ";";
							Records fastaLine_seqs = select_stmt(sql_command.c_str(), db);

							for (auto &fastaLine_seq : fastaLine_seqs) {
								auto fastaLine = fastaLine_seq[0];
								int align_start = fastaLine.find(kmer);
								int kmerend = align_start + kmer_size;
								string matched_seq_fasta = kmer;
								string matched_seq_fastq = kmer;
								int read_length = fastqLine.length();

								fastaLine = fastaLine.substr(align_start, fastaLine.length());
								// Is fastaLine shorter than read? then append next line of fasta so that length doesnt run out
								if (read_length < fastaLine.length()) {
									int next_fastaLine_id = stoi(fastaLine_id.at(0)) + 1;
									sql_command =
											"SELECT transcripts FROM ref_table WHERE id=" +
											to_string(next_fastaLine_id) + ";";
									Records next_fastaLine = select_stmt(sql_command.c_str(), db);
									string next_fastaLine_str = next_fastaLine.at(0).at(
											0); // at(0) twice because its a string inside a the Record vector, itself inside a Records vector
									fastaLine = fastaLine + next_fastaLine_str;
								}

								for (int i = kmer_size; i < read_length; i++) {

									int fastq_rel_pos = i;
									int fasta_rel_pos = fastq_rel_pos + align_start;

									int fasta_line_offset_test = is_offset(fastaLine, fasta_rel_pos, fastqLine,
									                                       fastq_rel_pos);
									int fastq_line_offset_test = is_offset(fastqLine, fastq_rel_pos, fastaLine,
									                                       fasta_rel_pos);

									if (fastqLine[fastq_rel_pos] == fastaLine[fasta_rel_pos]) {
										matched_seq_fastq += fastaLine[fastq_rel_pos];
										matched_seq_fasta += fastaLine[fastq_rel_pos];

									} else if (fasta_line_offset_test > 0) {
										int gapSize = fasta_line_offset_test; // returns size of gap
										auto x = fastaLine.substr(0, fasta_rel_pos);
										fastaLine =
												fastaLine.substr(0, fasta_rel_pos) +
												insert_char_x_times('-', gapSize) +
												fastaLine.substr(fasta_rel_pos, fastaLine.length());
										matched_seq_fasta += "-";
										matched_seq_fastq += fastqLine[fastq_rel_pos];

									} else if (fastq_line_offset_test > 0) {
										int gapSize = fastq_line_offset_test; // returns size of gap
										fastqLine =
												fastqLine.substr(0, fastq_rel_pos) +
												insert_char_x_times('-', gapSize) +
												fastqLine.substr(fastq_rel_pos, fastqLine.length());
										matched_seq_fastq += "-";
										matched_seq_fasta += fastaLine[fasta_rel_pos];

									} else { // mismatch
										matched_seq_fasta += fastaLine[fasta_rel_pos];
										matched_seq_fastq += fastqLine[fastq_rel_pos];
									}
								}


								double my_penalty = calculate_penalty(matched_seq_fasta, matched_seq_fastq);

								if (my_penalty < penalty_cutoff) {
									double my_similarity = calculate_similarity(matched_seq_fasta,
									                                            matched_seq_fastq);
									similarity_scores.push_back(my_similarity);
									if (verbosity_level > 2) {
										cout << "\nPenalty : " << to_string(my_penalty) << endl;
										cout << "Similarity : " << to_string(my_similarity) << endl;
										cout << matched_seq_fastq << endl;
										cout << matched_seq_fasta << endl;
									}

									string fastaLine_id_str = fastaLine_id.at(0);
									sql_command =
											"UPDATE ref_table SET quant = quant + 1 WHERE id=" + fastaLine_id_str +
											";";
									sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);
									aligned_counter += 1;
								}
							}
						}
					}
				}
			}
		}



		// Print quantification scores
		if (verbosity_level > 1) {
			sql_command = "SELECT header,SUM(quant) FROM ref_table GROUP BY header;";
			Records quant_scores = select_stmt(sql_command.c_str(), db);
			cout << "\nQuantification:" << endl;
			for (Record score : quant_scores) {
				cout << score.at(1) << "\t" << score.at(0) << endl;
			}
		}

		int similarity_scores_len = similarity_scores.size();
		double similarity_scores_sum;
		for (double score : similarity_scores) {
			similarity_scores_sum += score;
		}

		sql_command = "SELECT SUM(quant) FROM ref_table WHERE quant>0;";
		Records aligned_reads = select_stmt(sql_command.c_str(), db);
		double aligned_reads_rate = 100 * stod(aligned_reads.at(0).at(0)) / (double) read_counter;

		if (verbosity_level > 1) {
			cout << "\nAverage similarity score for the run: " << similarity_scores_sum / similarity_scores_len << "% "
			     << endl;
			cout << "Aligned Read rate: " << aligned_reads_rate << "\t nb aligned: " << aligned_counter << " / "
			     << read_counter << endl;
			cout << "% \t Time: ";
		}

		if (verbosity_level == 1) {
			cout << ::kmer_size << "," << ::max_gap_size << "," << ::penalty_cutoff << "," << ::mismatch_penalty << ","
			     << ::gap_opening_penalty << "," << ::gap_continuing_penalty << ","
			     << similarity_scores_sum / similarity_scores_len << "," << aligned_reads_rate << ",";

		}
	}
	sqlite3_close(db);

	auto end = chrono::steady_clock::now();
	auto time_diff = end - start;
	if (verbosity_level > 0) {
		cout << chrono::duration<double, milli>(time_diff).count() << endl;
	}

}



