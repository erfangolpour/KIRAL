#include <atomic>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <string>
#include <thread>

#include "align_thread.hpp"
#include "cli.hpp"
#include "helper.hpp"
#include "kir.hpp"
#include "types.hpp"

using namespace std;

int main(int argc, char const *argv[]) {
    if (argc < 2)
        return show_help(argv[0]);
    string command = argv[1];
    if (command == "-h" || command == "--help" || command == "help")
        return show_help(argv[0]);
    else if (command == "report") {
        if (argc < 3)
            return show_help(argv[0]);

        // Parse arguments
        string kir_id = "";
        string allele_id = "";
        int read_id = -1;
        int num_results = -1;
        string alignments_file = argv[2];
        for (int i = 3; i < argc; i++)
            if (string(argv[i]) == "--head")
                num_results = stoi(argv[++i]);
            else if (string(argv[i]) == "-r")
                read_id = stoi(argv[++i]);
            else if (string(argv[i]) == "-k")
                kir_id = argv[++i];
            else if (string(argv[i]) == "-a")
                allele_id = argv[++i];
        show_report(alignments_file, kir_id, allele_id, read_id, num_results);
    } else if (command == "align") {
        if (argc < 4)
            return show_help(argv[0]);

        // Parse arguments
        string kirs_file = argv[2];
        string reads_file = argv[3];
        string method = "regional";
        int num_representatives = 1;
        bool inc_pair = false;
        int n_threads = thread::hardware_concurrency();
        string output_file = "";
        for (int i = 4; i < argc; i++)
            if (string(argv[i]) == "--method") {
                method = argv[++i];
                if (method != "naive" && method != "regional" && method != "categorical")
                    return show_help(argv[0]);
            } else if (string(argv[i]) == "-r")
                num_representatives = stoi(argv[++i]);
            else if (string(argv[i]) == "--pair")
                inc_pair = true;
            else if (string(argv[i]) == "-t")
                n_threads = stoi(argv[++i]);
            else if (string(argv[i]) == "-o")
                output_file = argv[++i];
        cout << "[+] Using " << n_threads << " thread(s)." << endl;

        // Load data
        unordered_map<string, unordered_map<string, string>> kirs = load_kirs(kirs_file);
        unordered_map<int, string> reads = load_reads(reads_file);
        cout << "[+] Loaded " << reads.size() << " reads." << endl;

        // Extract reads to fasta file with new indexing
        cout << "[*] Extracting reads to fasta file with new indexing..." << flush;
        string reads_fasta_file = reindex_reads(reads);
        cout << "\r[✓]" << endl;

        // Mutex for thread-safe access to all_alignments
        unordered_map<string, unordered_map<string, vector<ReadAlignment>>> all_alignments;
        vector<thread> threads;
        mutex mtx;

        // Atomic variable to track progress
        atomic<int> progress(0);
        int total_genes = kirs.size();
        int bar_width = 70;

        // Function to display progress bar
        auto display_progress = [&]() {
            while (progress < total_genes) {
                float progress_ratio = static_cast<float>(progress) / total_genes;
                int pos = bar_width * progress_ratio;
                cout << "\r[";
                for (int i = 0; i < bar_width; ++i)
                    cout << (i < pos ? "=" : (i == pos ? ">" : " "));
                cout << "] " << int(progress_ratio * 100.0) << " %";
                cout.flush();
                this_thread::sleep_for(chrono::milliseconds(100));
            }
            cout << "\r[";
            for (int i = 0; i < bar_width; ++i)
                cout << "=";
            cout << "] 100 %\n";
        };

        // Perform alignment
        if (method == "naive") {
            auto gene_it = kirs.begin();
            mutex gene_it_mtx;
            cout << "[*] Performing naive alignment..." << endl;
            thread progress_thread(display_progress);

            for (int i = 0; i < n_threads; ++i)
                threads.push_back(thread(naive_align, i, ref(kirs), ref(reads_fasta_file), ref(all_alignments), ref(mtx), ref(gene_it), ref(gene_it_mtx), ref(progress), n_threads));

            // Wait for all threads to finish
            for (auto &t : threads)
                t.join();

            // Wait for progress thread to finish
            progress_thread.join();
        } else {
            // Regional and categorical alignment both require a first pass to extract representative alleles
            cout << "[*] Extracting " << num_representatives << " representative allele(s) per gene..." << flush;
            string representatives_file = extract_representatives(kirs, num_representatives);
            cout << "\r[✓]" << endl;

            cout << "[*] Performing initial alignment with representative alleles..." << flush;
            auto first_pass_results = align_minimap(representatives_file, reads_fasta_file, n_threads);
            cleanup(representatives_file);
            cout << "\r[✓]" << endl;

            auto gene_it = first_pass_results.begin();
            mutex gene_it_mtx;

            cout << "[*] Performing " << method << " alignment on " << total_genes << " gene(s)..." << endl;
            thread progress_thread(display_progress);

            auto method_func = "regional" ? regional_align : categorical_align;
            for (int i = 0; i < n_threads; ++i)
                threads.push_back(thread(method_func, i, ref(kirs), ref(reads), ref(first_pass_results), ref(all_alignments), ref(mtx), ref(gene_it), ref(gene_it_mtx), ref(progress), inc_pair, n_threads));

            // Wait for all threads to finish
            for (auto &t : threads)
                t.join();

            // Wait for progress thread to finish
            progress_thread.join();
        }

        // Cleanup intermediate files
        cleanup(reads_fasta_file);

        // Sort result lines by adding them to a heap first
        // uint total_matches = 0;
        // priority_queue<string, vector<string>, greater<string>> results;  // min heap
        // for (const auto &gene : all_alignments) {
        //     uint total = 0;
        //     for (const auto &allele : gene.second)
        //         total += allele.second.size();
        //     total_matches += total;
        //     results.push(gene.first + ":" + string(10 - gene.first.size(), ' ') + to_string(total) + " matches found in " + to_string(gene.second.size()) + " alleles");
        // }

        // // Report results
        // cout << "\n[+] Results:" << endl;
        // while (!results.empty()) {
        //     cout << results.top() << endl;
        //     results.pop();
        // }
        // cout << "[+] Total matches: " << total_matches << endl;

        if (!output_file.empty()) {
            ofstream out_file(output_file);
            if (out_file.is_open()) {
                for (const auto &gene : all_alignments)
                    for (const auto &allele : gene.second)
                        for (const auto &alignment : allele.second)
                            out_file << alignment.read_id << "\t" << alignment.kir_id << "\t" << alignment.allele_id << "\t" << (alignment.reversed ? "1" : "0") << "\t" << alignment.cost << "\t" << alignment.read_start << "\t" << alignment.read_end << "\t" << alignment.query_start << "\t" << alignment.query_end << "\t" << alignment.cigar << endl;
                out_file.close();
            } else
                cerr << "[-] Error: Unable to open file " << output_file << " for writing." << endl;
        }
        cout << "[+] Results saved to " << output_file << endl;

    } else
        return show_help(argv[0]);

    // Return success
    return 0;
}
