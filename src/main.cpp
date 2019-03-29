#include <chrono>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <sqlite3.h>
#include "cxxopts.hpp"

using namespace std;
auto start = chrono::steady_clock::now();

string fasta_path = "example_data/mini_cdna_CRCh37";
string fastq_path = "example_data/sim_1.fq";
//string fastq_path = "example_data/small/myfastq.fastq";
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

using Record = vector<string>;
using Records = vector<Record>;

int select_callback(void *p_data, int num_fields, char **p_fields, char **p_col_names) {
	Records *records = static_cast<Records *>(p_data);
	try {
		records->emplace_back(p_fields, p_fields + num_fields);
	}
	catch (...) {
		// abort select on failure, don't let exception propogate thru sqlite3 call-stack
		return 1;
	}
	return 0;
}

Records select_stmt(const char *stmt, sqlite3 *db) {
	Records records;
	char *errmsg;
	int ret = sqlite3_exec(db, stmt, select_callback, &records, &errmsg);
	if (ret != SQLITE_OK) {
		std::cerr << "Error in select statement " << stmt << "[" << errmsg << "]\n";
	}

	return records;
}

string insert_char_x_times(char input_char, int times) {
	string new_string;
	for (int i = 1; i <= times; ++i) {
		new_string += input_char;
	}
	return new_string;
}

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


static int callback(void *NotUsed, int argc, char **argv, char **azColName) {
	int i;
	for (i = 0; i < argc; i++) {
		cout << azColName[i] << " = " << (argv[i] ? argv[i] : "NULL") << "\n";
	}
	return 0;
}


int main(int argc, char *argv[]) {
	cxxopts::Options options("MyProgram", "One line description of MyProgram");
	options.add_options()
			("skip-db", "Skip fasta parsing and database construction", cxxopts::value<bool>()->default_value("false"))
			("skip-quant", "Skip quant, just build index database", cxxopts::value<bool>()->default_value("false"))
			("v,verbose",
			 "Set verbosity level, 0 is no output, 1 is CSV of results, 2 is quant scores, 3 is read alignment",
			 cxxopts::value<int>()->default_value("3"))
			("k,kmer-size", "kmer size to search fasta", cxxopts::value<int>())
			("c,cutoff", "cutoff above which read with this penalty won't be quantified", cxxopts::value<int>())
			("g,max-gap", "maximum gapsize to look for gaps", cxxopts::value<int>())
			("m,mismatch", "penalty for nucleotide mismatch", cxxopts::value<int>())
			("o,gap-open", "penalty for opening a gap", cxxopts::value<int>())
			("e,gap-extend", "penalty for extending a gap", cxxopts::value<int>());

	auto result = options.parse(argc, argv);
	::skip_quant = result["skip-quant"].as<bool>();
	::skip_database_creation = result["skip-db"].as<bool>();
	::verbosity_level = result["verbose"].as<int>();
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
		sql_command = "CREATE TABLE IF NOT EXISTS ref_table(id INTEGER PRIMARY KEY  AUTOINCREMENT, quant INTEGER, header TEXT, transcripts VARCHAR(80));";
		sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);

		sqlite3_exec(db, "PRAGMA cache_size=10000", NULL, NULL, NULL);
		sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);


		string line;
		string current_header;
		while (getline(fasta_file, line)) {
			if (line[0] == '>') {
				line = line.substr(1, line.size());
				line = regex_replace(line, regex("^ "), "");
				line = regex_replace(line, regex(" +"), " ");
				current_header = line;

			} else {
				sql_command =
						"INSERT INTO ref_table (transcripts, quant, header) VALUES ('" + line + "',0,'" +
						current_header +
						"');";
				sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);
			}
		}
		sqlite3_exec(db, "END TRANSACTION", NULL, NULL, NULL);

		// Index database for log n lookup times
		sql_command = "CREATE INDEX transcript_index ON ref_table (transcripts);CREATE INDEX id_index ON ref_table (id);";
		sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);
	}

	if (::skip_quant == false) {
		int line_counter = 0;
		int read_counter = 0;
		string fastqLine;
		while (getline(fastq_file, fastqLine)) {
			line_counter += 1;
			if (line_counter % 4 == 2) {
				read_counter += 1;
				string kmer;
				kmer = fastqLine.substr(0, kmer_size);

				// Get id of database row where kmer exists (less expensive option because O(log n) )
				sql_command = "SELECT id FROM ref_table WHERE transcripts LIKE '" + kmer + "%';";
				Records fastaLine_ids = select_stmt(sql_command.c_str(), db);

			if (fastaLine_ids.data() == 0x0) {
				// using more expensive O(n) kmer search if can't find kmer with cheaper option
				sql_command = "SELECT id FROM ref_table WHERE transcripts LIKE '%" + kmer + "%';";
				fastaLine_ids = select_stmt(sql_command.c_str(), db);
			}
			for (auto fastaLine_id : fastaLine_ids) {
				sql_command = "SELECT transcripts FROM ref_table WHERE id=" + (string) fastaLine_id[0] + ";";
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
									"SELECT transcripts FROM ref_table WHERE id=" + to_string(next_fastaLine_id) + ";";
							Records next_fastaLine = select_stmt(sql_command.c_str(), db);
							string next_fastaLine_str = next_fastaLine.at(0).at(
									0); // at(0) twice because its a string inside a the Record vector, itself inside a Records vector
							fastaLine = fastaLine + next_fastaLine_str;
						}

						for (int i = kmer_size; i < read_length; i++) {

							int fastq_rel_pos = i;
							int fasta_rel_pos = fastq_rel_pos + align_start;

							int fasta_line_offset_test = is_offset(fastaLine, fasta_rel_pos, fastqLine, fastq_rel_pos);
							int fastq_line_offset_test = is_offset(fastqLine, fastq_rel_pos, fastaLine, fasta_rel_pos);

							if (fastqLine[fastq_rel_pos] == fastaLine[fasta_rel_pos]) {
								matched_seq_fastq += fastaLine[fastq_rel_pos];
								matched_seq_fasta += fastaLine[fastq_rel_pos];

							} else if (fasta_line_offset_test > 0) {
								int gapSize = fasta_line_offset_test; // returns size of gap
								auto x = fastaLine.substr(0, fasta_rel_pos);
								fastaLine = fastaLine.substr(0, fasta_rel_pos) + insert_char_x_times('-', gapSize) +
								            fastaLine.substr(fasta_rel_pos, fastaLine.length());
								matched_seq_fasta += "-";
								matched_seq_fastq += fastqLine[fastq_rel_pos];

							} else if (fastq_line_offset_test > 0) {
								int gapSize = fastq_line_offset_test; // returns size of gap
								fastqLine = fastqLine.substr(0, fastq_rel_pos) + insert_char_x_times('-', gapSize) +
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
							double my_similarity = calculate_similarity(matched_seq_fasta, matched_seq_fastq);
							similarity_scores.push_back(my_similarity);
							if (verbosity_level > 2) {
								cout << "\nPenalty : " << to_string(my_penalty) << endl;
								cout << "Similarity : " << to_string(my_similarity) << endl;
								cout << matched_seq_fastq << endl;
								cout << matched_seq_fasta << endl;
							}

							string fastaLine_id_str = fastaLine_id.at(0);
							sql_command = "UPDATE ref_table SET quant = quant + 1 WHERE id=" + fastaLine_id_str + ";";
							sqlite3_exec(db, sql_command.c_str(), callback, 0, NULL);
						}
					}
				}
			}
		}


		// Print quantification scores
		if (verbosity_level > 1) {
			sql_command = "SELECT header,SUM(quant) FROM ref_table GROUP BY header;";
			Records quant_scores = select_stmt(sql_command.c_str(), db);
			cout << "\nScores:" << endl;
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
			     << similarity_scores_sum << "/" << similarity_scores_len << endl;
			cout << "Aligned Read rate: " << aligned_reads_rate << "% \t Time: ";
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



