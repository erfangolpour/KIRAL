#ifndef ALIGN_THREAD_H
#define ALIGN_THREAD_H

#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <set>
#include <vector>
#include <mutex>
#include <atomic>
#include <algorithm>

#include "helper.hpp"
#include "types.hpp"
#include "kir.hpp"

using namespace std;

void naive_align(int thread_id, unordered_map<string, unordered_map<string, string>> &kirs, const string &reads_file, unordered_map<string, unordered_map<string, vector<ReadAlignment>>> &all_alignments, mutex &mtx, unordered_map<string, unordered_map<string, string>>::iterator &gene_it, mutex &gene_it_mtx, atomic<int> &progress, int n_threads)
{
    string alleles_fasta_file = "alleles_" + to_string(thread_id) + ".fa";

    while (true)
    {
        // get slice of [i * len(kirs) / n_threads, (i + 1) * len(kirs) / n_threads]
        unordered_map<string, unordered_map<string, string>> thread_kirs;
        {
            lock_guard<mutex> lock(gene_it_mtx);
            if (gene_it == kirs.end())
                break;
            for (int i = 0; i < (int)kirs.size() / n_threads; ++i)
            {
                if (gene_it == kirs.end())
                    break;
                thread_kirs[gene_it->first] = gene_it->second;
                ++gene_it;
            }
        }

        ofstream alleles_fasta(alleles_fasta_file);
        for (const auto &gene : thread_kirs)
            for (const auto &allele : gene.second)
                alleles_fasta << ">" << gene.first << "." << allele.first << "\n"
                              << allele.second << "\n";
        alleles_fasta.close();

        auto results = align_minimap(alleles_fasta_file, reads_file, n_threads);

        {
            lock_guard<mutex> lock(mtx);
            for (const auto &gene_alignments : results)
                for (const auto &allele_alignments : gene_alignments.second)
                    all_alignments[gene_alignments.first][allele_alignments.first].insert(
                        all_alignments[gene_alignments.first][allele_alignments.first].end(),
                        allele_alignments.second.begin(), allele_alignments.second.end());
        }

        // Update progress
        progress += thread_kirs.size();
    }

    // Cleanup intermediate files
    cleanup(alleles_fasta_file);
}

void regional_align(int thread_id, unordered_map<string, unordered_map<string, string>> &kirs, unordered_map<int, string> &reads, unordered_map<string, unordered_map<string, vector<ReadAlignment>>> &first_pass_results, unordered_map<string, unordered_map<string, vector<ReadAlignment>>> &all_alignments, mutex &mtx, unordered_map<string, unordered_map<string, vector<ReadAlignment>>>::iterator &gene_it, mutex &gene_it_mtx, atomic<int> &progress, bool inc_pair, int n_threads)
{
    string reads_fasta_file = "gene_reads_" + to_string(thread_id) + ".fa";
    string alleles_fasta_file = "alleles_" + to_string(thread_id) + ".fa";

    while (true)
    {
        // Get the next gene to process
        string gene_name;
        {
            lock_guard<mutex> lock(gene_it_mtx);
            if (gene_it == first_pass_results.end())
                break;
            gene_name = gene_it->first;
            ++gene_it;
        }

        for (const auto &allele_alignments : first_pass_results[gene_name])
        {
            vector<Region> regions;
            // TODO: Make this a parameter, for now
            int region_buffer = reads.begin()->second.size() * 100; // Assume all reads have the same length

            // Find reads that belong to this region
            for (const auto &alignment : allele_alignments.second)
            {
                Region region(alignment.query_start, alignment.query_end, region_buffer);
                region.add_read(alignment.read_id);
                if (inc_pair)
                    region.add_read(get_pair_id(alignment.read_id));

                auto it = find(regions.begin(), regions.end(), region); // two overlapping regions are considered equal
                if (it != regions.end())
                    it->merge(region);
                else
                    regions.push_back(region);
            }

            // Filter out regions that are too small and add their reads to a single region
            Region common_region(numeric_limits<int>::max(), 0, region_buffer);
            for (auto it = regions.begin(); it != regions.end();)
            {
                if (it->reads.size() < 3) // TODO: Make this a parameter
                {
                    common_region.merge(*it);
                    it = regions.erase(it);
                }
                else
                    it++;
            }
            regions.push_back(common_region);

            for (const auto &region : regions)
            {
                if (region.reads.empty())
                    continue;
    
                // Extract reads that belong to this region
                ofstream gene_reads(reads_fasta_file);
                for (const auto &read_id : region.reads)
                    gene_reads << ">" << read_id << "\n"
                               << reads[read_id] << "\n";
                gene_reads.close();

                // generate a temporary allele fasta where all alleles are trimmed to the region
                ofstream alleles_fasta(alleles_fasta_file);
                for (const auto &allele : kirs[gene_name])
                {
                    string trimmed_allele = allele.second.substr(min(region.start, (int)allele.second.size() - 1), region.end - region.start);
                    alleles_fasta << ">" << gene_name << "." << allele.first << "\n"
                                  << trimmed_allele << "\n";
                }
                alleles_fasta.close();

                auto second_pass_results = align_minimap(alleles_fasta_file, reads_fasta_file, n_threads);
                if (second_pass_results.find(gene_name) == second_pass_results.end())
                    continue;

                // Merge the results with the global alignments map
                for (auto &allele_alignments : second_pass_results[gene_name])
                    for (auto &alignment : allele_alignments.second)
                    {
                        alignment.query_start += region.start;
                        alignment.query_end = min(alignment.query_end + region.start, (int)kirs[gene_name][allele_alignments.first].size());
                        all_alignments[gene_name][allele_alignments.first].push_back(alignment);
                    }
            }
        }

        // Update progress
        progress++;
    }

    // Cleanup intermediate files
    cleanup(reads_fasta_file);
    cleanup(alleles_fasta_file);
}

void categorical_align(int thread_id, unordered_map<string, unordered_map<string, string>> &kirs, unordered_map<int, string> &reads, unordered_map<string, unordered_map<string, vector<ReadAlignment>>> &first_pass_results, unordered_map<string, unordered_map<string, vector<ReadAlignment>>> &all_alignments, mutex &mtx, unordered_map<string, unordered_map<string, vector<ReadAlignment>>>::iterator &gene_it, mutex &gene_it_mtx, atomic<int> &progress, bool inc_pair, int n_threads)
{
    string reads_fasta_file = "reads_" + to_string(thread_id) + ".fa";
    string alleles_fasta_file = "alleles_" + to_string(thread_id) + ".fa";

    while (true)
    {
        // Get the next gene to process
        string gene_name;
        {
            lock_guard<mutex> lock(gene_it_mtx);
            if (gene_it == first_pass_results.end())
                break;
            gene_name = gene_it->first;
            ++gene_it;
        }

        // extract the ID of the reads that aligned to this gene
        set<int> read_ids;
        for (const auto &allele_alignments : first_pass_results[gene_name])
        {
            for (const auto &match : allele_alignments.second)
            {
                read_ids.insert(match.read_id);
                if (inc_pair)
                    read_ids.insert(get_pair_id(match.read_id));
            }
        }

        ofstream gene_reads(reads_fasta_file);
        for (const auto &id : read_ids)
            gene_reads << ">" << id << "\n"
                       << reads[id] << "\n";
        gene_reads.close();

        ofstream alleles_fasta(alleles_fasta_file);
        for (const auto &allele : kirs[gene_name])
            alleles_fasta << ">" << gene_name << "." << allele.first << "\n"
                          << allele.second << "\n";
        alleles_fasta.close();

        auto second_pass_results = align_minimap(alleles_fasta_file, reads_fasta_file, n_threads);
        if (second_pass_results.find(gene_name) == second_pass_results.end())
            // No matches in the entire gene, we should ideally never reach this point, or something must have gone wrong
            continue;

        {
            lock_guard<mutex> lock(mtx);
            for (const auto &allele_alignments : second_pass_results[gene_name])
                all_alignments[gene_name][allele_alignments.first].insert(
                    all_alignments[gene_name][allele_alignments.first].end(),
                    allele_alignments.second.begin(), allele_alignments.second.end());
        }

        // Update progress
        progress++;
    }

    // Cleanup intermediate files
    cleanup(reads_fasta_file);
    cleanup(alleles_fasta_file);
}

#endif