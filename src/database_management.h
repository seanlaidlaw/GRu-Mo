//
// Created by Sean Laidlaw on 2019-04-03.
//

#ifndef GRU_MO_DATABASE_MANAGEMENT_H
#define GRU_MO_DATABASE_MANAGEMENT_H

using namespace std;
using Record = vector<string>;
using Records = vector<Record>;

int select_callback(void *p_data, int num_fields, char **p_fields, char **p_col_names) {
	Records *records = static_cast<Records *>(p_data);
	try {
		records->emplace_back(p_fields, p_fields + num_fields);
	}
	catch (...) {
		// abort select on failure
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

static int callback(void *NotUsed, int argc, char **argv, char **azColName) {
	int i;
	for (i = 0; i < argc; i++) {
		cout << azColName[i] << " = " << (argv[i] ? argv[i] : "NULL") << "\n";
	}
	return 0;
}

#endif //GRU_MO_DATABASE_MANAGEMENT_H
