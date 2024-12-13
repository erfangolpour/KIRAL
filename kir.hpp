#ifndef KIR_H
#define KIR_H

#include <vector>
#include <iostream>
#include <zlib.h>
#include <fstream>
#include <string>
#include <unordered_map>

#include "types.hpp"
#include "helper.hpp"
#include "minimap2/kseq.h"
#include "minimap2/minimap.h"

using namespace std;

KSEQ_INIT(gzFile, gzread)

int get_pair_id(int read_id)
{
    return read_id - 1 + !(read_id % 2) * 2;
}

unordered_map<string, unordered_map<string, string>> load_kirs(const string &database)
{
    unordered_map<string, unordered_map<string, string>> kirs;
    gzFile kirFileIn = expect(gzopen(database.c_str(), "r"), "Failed to open database file");
    kseq_t *seq = kseq_init(kirFileIn);
    while (kseq_read(seq) >= 0)
    {
        string kir_id = seq->name.s;
        auto pos = kir_id.find('.');
        string gene_id = kir_id.substr(0, pos);
        string allele_id = kir_id.substr(pos + 1);
        kirs[gene_id][allele_id] = seq->seq.s;
    }
    kseq_destroy(seq);
    gzclose(kirFileIn);
    return kirs;
}

unordered_map<int, string> load_reads(const string &reads)
{
    int r_id = 0;
    unordered_map<int, string> reads_map;
    gzFile readsFileIn = expect(gzopen(reads.c_str(), "r"), "Failed to open reads file");
    kseq_t *seq = kseq_init(readsFileIn);
    while (kseq_read(seq) >= 0)
        reads_map[r_id++] = seq->seq.s;
    kseq_destroy(seq);
    gzclose(readsFileIn);
    return reads_map;
}

const string reindex_reads(const unordered_map<int, string> &reads_map)
{
    string reads_fasta = "reads_indexed.fa";
    ofstream reads_fasta_file(reads_fasta);
    for (const auto &read : reads_map)
        reads_fasta_file << ">" << read.first << "\n"
                         << read.second << "\n";
    reads_fasta_file.close();
    return reads_fasta;
}

const string extract_representatives(const unordered_map<string, unordered_map<string, string>> &kirs, int num_representatives)
{
    string representatives_file = "representatives.fa";
    ofstream representatives(representatives_file);
    for (const auto &gene : kirs)
        // Choose n random alleles to represent the gene
        for (int i = 0; i < num_representatives; i++)
        {
            auto allele = next(gene.second.begin(), rand() % gene.second.size());
            representatives << ">" << gene.first << "." << allele->first << "\n"
                            << allele->second << "\n";
        }
    representatives.close();
    return representatives_file;
}

unordered_map<string, unordered_map<string, vector<ReadAlignment>>> align_minimap(
    const string &kirdb, const string &reads, int n_threads, int max_num_mismatches = 5)
{
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    expect(!mm_set_opt(0, &iopt, &mopt), "Failed to set minimap options");
    mopt.flag |= MM_F_SR;         // -x sr: Short-read preset
    mopt.flag |= MM_F_CIGAR;      // -c: Calculate CIGAR strings
    mopt.flag |= MM_F_ALL_CHAINS; // -P: Retain all chains and don’t attempt to set primary chains

    // Open query file
    gzFile readsFile = expect(gzopen(reads.c_str(), "r"), "Failed to open reads file");
    kseq_t *ks = kseq_init(readsFile);

    unordered_map<string, unordered_map<string, vector<ReadAlignment>>> alignments;
    mm_idx_reader_t *r = mm_idx_reader_open(kirdb.c_str(), &iopt, NULL);
    mm_idx_t *mi; // minimap2 index
    while ((mi = mm_idx_reader_read(r, n_threads)) != 0)
    {
        mm_mapopt_update(&mopt, mi);      // this sets the maximum minimizer occurrence
        mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
        gzrewind(readsFile);
        kseq_rewind(ks);
        while (kseq_read(ks) >= 0)
        {
            mm_reg1_t *reg;
            int j, n_reg;
            reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, NULL); // get all hits for the query
            for (j = 0; j < n_reg; ++j)
            {
                int num_mismatches = reg[j].blen - reg[j].mlen + reg[j].p->n_ambi;
                if (num_mismatches + (int)ks->seq.l - (reg[j].re - reg[j].rs) > max_num_mismatches)
                    continue;
                string kir_key = string(mi->seq[reg[j].rid].name);
                auto pos = kir_key.find('.');
                string gene_key = kir_key.substr(0, pos);
                string allele_key = kir_key.substr(pos + 1);
                string cigar;
                for (uint32_t k = 0; k < reg[j].p->n_cigar; k++) // this gives the CIGAR in the aligned regions. NO soft/hard clippings!
                    cigar += to_string(reg[j].p->cigar[k] >> 4) + MM_CIGAR_STR[reg[j].p->cigar[k] & 0xf];
                ReadAlignment match = {stoi(ks->name.s), gene_key, allele_key, reg[j].rev != 0, num_mismatches, reg[j].rs, reg[j].re, reg[j].qs, reg[j].qe, cigar};
                alignments[gene_key][allele_key].push_back(match);
                free(reg[j].p);
            }
            free(reg);
        }
        mm_tbuf_destroy(tbuf); // deallocate the thread buffer
        mm_idx_destroy(mi);    // deallocate the index
    }

    mm_idx_reader_close(r); // close the index reader
    kseq_destroy(ks);       // close the query file
    gzclose(readsFile);     // close the query file
    
    return alignments;
}

#endif