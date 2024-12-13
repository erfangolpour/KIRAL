#ifndef CLI_H
#define CLI_H

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

#include "types.hpp"

using namespace std;

void print_alignment(const ReadAlignment &alignment)
{
    cout << "Read ID: " << alignment.read_id << endl;
    cout << "KIR ID: " << alignment.kir_id << endl;
    cout << "Allele ID: " << alignment.allele_id << endl;
    cout << "Reversed: " << (alignment.reversed ? "Yes" : "No") << endl;
    cout << "Cost: " << alignment.cost << endl;
    cout << "Read Start: " << alignment.read_start << endl;
    cout << "Read End: " << alignment.read_end << endl;
    cout << "Query Start: " << alignment.query_start << endl;
    cout << "Query End: " << alignment.query_end << endl;
    cout << "CIGAR: " << alignment.cigar << endl;
    cout << endl;
}

void show_report(const string &alignments_file, const string &kir_id, const string &allele_id, int read_id, int num_results)
{
    ifstream alignments(alignments_file);
    if (!alignments)
    {
        cerr << "[x] Failed to open alignments file " << alignments_file << endl;
        return;
    }
    string line;
    while (num_results-- && getline(alignments, line))
    {
        ReadAlignment alignment;
        // split the line by tabs
        size_t pos = 0;
        vector<string> fields;
        while ((pos = line.find('\t')) != string::npos)
        {
            fields.push_back(line.substr(0, pos));
            line.erase(0, pos + 1);
        }
        fields.push_back(line);
        // parse the fields
        alignment.read_id = stoi(fields[0]);
        if (read_id != -1 && alignment.read_id != read_id)
            continue;
        alignment.kir_id = fields[1];
        if (!kir_id.empty() && alignment.kir_id != kir_id)
            continue;
        alignment.allele_id = fields[2];
        if (!allele_id.empty() && alignment.allele_id != allele_id)
            continue;
        alignment.reversed = fields[3] == "1";
        alignment.cost = stoi(fields[4]);
        alignment.read_start = stoi(fields[5]);
        alignment.read_end = stoi(fields[6]);
        alignment.query_start = stoi(fields[7]);
        alignment.query_end = stoi(fields[8]);
        alignment.cigar = fields[9];
        print_alignment(alignment);
    }
}

/* Function to print the help message */
int show_help(const string &program)
{
    cerr << "Usage: " << program << " <command> [options]\n"
         << endl;
    cerr << "Commands:" << endl;

    cerr << "\talign <database> <reads> [--method <method_name>] [-r <num_representatives>] [--pair] [-t <threads>] [-o <output_file>]" << endl;
    cerr << "\t\tAligns reads to the database and reports the results." << endl;
    cerr << "\tOptions:" << endl;
    cerr << "\t\t--method <method_name>\n"
         << "\t\t\tAlignment method to use. Options are `naive`, `regional`, and `categorical`. Default is `regional`." << endl;
    cerr << "\t\t-r <num_representatives>\n"
         << "\t\t\tNumber of representative alleles per gene used in `regional` and `categorical` alignment. Default is 1." << endl;
    cerr << "\t\t--pair\n"
         << "\t\t\tWhen performing `regional` or `categorical` alignment, for each read aligned in the first pass, also include its pair in the second pass." << endl;
    cerr << "\t\t-t <threads>\n"
         << "\t\t\tNumber of threads to use. Default is the number of hardware threads." << endl;
    cerr << "\t\t-o <output_file>\n"
         << "\t\t\tOutput file to write the results to." << endl;

    cerr << "\n\treport <alignments_file> [--head <num_results>] [-r <read read_id>] [-k <KIR read_id>] [-a <allele read_id>]" << endl;
    cerr << "\t\tReports the results from a previously generated alignments file." << endl;
    cerr << "\tOptions:" << endl;
    cerr << "\t\t--head <num_results>\n"
         << "\t\t\tShow only the first <num_results> results." << endl;
    cerr << "\t\t-r <read_id>\n"
         << "\t\t\tShow only alignments for the specified read." << endl;
    cerr << "\t\t-k <KIR read_id>\n"
         << "\t\t\tShow only alignments for the specified KIR." << endl;
    cerr << "\t\t-a <allele read_id>\n"
         << "\t\t\tShow only alignments for the specified allele." << endl;
    return 1; // Return error code
}

#endif