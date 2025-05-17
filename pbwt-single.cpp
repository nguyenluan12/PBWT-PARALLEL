#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <fstream>
#include <iterator>
#include <bits/stdc++.h>
#include "util.h"

using namespace std;
using namespace chrono;

void build_prefix_and_divergence_arrays(const vector<vector<int> > &X, vector<int> &ppa, vector<int> &div) {
    // Number of haplotypes
    int M = X.size();

    // Number of variable sites
    int N = X[0].size();

    // Initialize positional prefix array and divergence array
    ppa.resize(M);
    div.resize(M, 0);

    // Fill initial ppa with indices
    for (int i = 0; i < M; ++i) {
        ppa[i] = i;
    }

    // Iterate over variants
    for (int k = 0; k < N; ++k) {
        // Temporary vectors to store intermediates
        vector<int> a, b, d, e;
        int p = k + 1, q = k + 1;

        // Iterate over haplotypes in reverse prefix sorted order
        for (int i = 0; i < M; ++i) {
            int index = ppa[i];
            int match_start = div[i];
            const vector<int> &haplotype = X[index];
            int allele = haplotype[k];

            // Update intermediates
            if (match_start > p) {
                p = match_start;
            }
            if (match_start > q) {
                q = match_start;
            }

            if (allele == 0) {
                a.push_back(index);
                d.push_back(p);
                p = 0;
            } else {
                b.push_back(index);
                e.push_back(q);
                q = 0;
            }
        }

        // Construct the new arrays for k+1 by concatenating intermediates
        ppa.clear();
        ppa.reserve(a.size() + b.size());
        ppa.insert(ppa.end(), a.begin(), a.end());
        ppa.insert(ppa.end(), b.begin(), b.end());

        div.clear();
        div.reserve(d.size() + e.size());
        div.insert(div.end(), d.begin(), d.end());
        div.insert(div.end(), e.begin(), e.end());
    }
}

vector<vector<int> > report_long_matches(const vector<vector<int> > &X, int L) {
    int M = X.size(); // Number of haplotypes
    int N = X[0].size(); // Number of variable sites

    // Initialize positional prefix array and divergence array
    vector<int> ppa(M);
    iota(ppa.begin(), ppa.end(), 0); // Fill ppa with 0, 1, 2, ..., M-1
    vector<int> div(M, 0);

    vector<vector<int> > result; // To store results

    // Iterate over variants
    for (int k = 0; k < N; ++k) {
        vector<int> a, b, d, e;
        int p = k + 1;
        int q = k + 1;
        vector<int> ma, mb;

        // Iterate over haplotypes in reverse prefix sorted order
        for (size_t i = 0; i < ppa.size(); ++i) {
            int index = ppa[i];
            int match_start = div[i];

            // Report matches
            if (match_start > k - L) {
                if (!ma.empty() && !mb.empty()) {
                    // Store position k, ma, and mb in the result vector
                    vector<int> result_entry;
                    result_entry.push_back(k);
                    result_entry.insert(result_entry.end(), ma.begin(), ma.end());
                    result_entry.push_back(-1); // Separator for mb
                    result_entry.insert(result_entry.end(), mb.begin(), mb.end());
                    result.push_back(result_entry);
                }
                ma.clear();
                mb.clear();
            }

            // Current haplotype
            const vector<int> &haplotype = X[index];
            int allele = haplotype[k];

            // Update intermediates
            if (match_start > p) p = match_start;
            if (match_start > q) q = match_start;

            // Update intermediates
            if (allele == 0) {
                a.push_back(index);
                d.push_back(p);
                p = 0;
                ma.push_back(index);
            } else {
                b.push_back(index);
                e.push_back(q);
                q = 0;
                mb.push_back(index);
            }
        }

        // Report any remaining matches including final haplotype
        if (!ma.empty() && !mb.empty()) {
            vector<int> result_entry;
            result_entry.push_back(k);
            result_entry.insert(result_entry.end(), ma.begin(), ma.end());
            result_entry.push_back(-1); // Separator for mb
            result_entry.insert(result_entry.end(), mb.begin(), mb.end());
            result.push_back(result_entry);
        }

        // Construct the new arrays for k+1
        if (k < N - 1) {
            vector<int> new_ppa;
            vector<int> new_div;

            new_ppa.reserve(a.size() + b.size());
            new_div.reserve(d.size() + e.size());

            new_ppa.insert(new_ppa.end(), a.begin(), a.end());
            new_ppa.insert(new_ppa.end(), b.begin(), b.end());

            new_div.insert(new_div.end(), d.begin(), d.end());
            new_div.insert(new_div.end(), e.begin(), e.end());

            ppa = move(new_ppa);
            div = move(new_div);
        }
    }

    return result;
}

int main(int argc, char *argv[]) {
    string f = argv[1];
    printf("Read file: %s\n", f.c_str());
    vector<vector<int> > X = read_hap(f);

    const int retry = 1;
    const int L = 4;
    printf("Test size: %lu haps x %lu sites\n", X.size(), X[0].size());
    printf("Retry: %d\n", retry);
    printf("Matched length: %d\n", L);

    vector<int> ppa;
    vector<int> div;
    signed long int duration_total_build = 0;
    high_resolution_clock::time_point start, stop;
    vector<vector<vector<int> > > res;
    for (int j = 0; j < retry; ++j) {
        start = high_resolution_clock::now();
        build_prefix_and_divergence_arrays(X, ppa, div);
        stop = high_resolution_clock::now();
        duration_total_build += duration_cast<microseconds>(stop - start).count();
    }
    printf("1 thread(s), build time: %ld us\n", duration_total_build / retry);

    signed long int duration_total_match = 0;
    for (int j = 0; j < retry; ++j) {
        start = high_resolution_clock::now();
        vector<vector<int> > matches = report_long_matches(X, L);
        stop = high_resolution_clock::now();
        duration_total_match += duration_cast<microseconds>(stop - start).count();
        print_matches(matches);
    }
    printf("1 thread(s), match time: %ld us\n", duration_total_match / retry);
    printf("1 thread(s), total time: %ld us\n\n", duration_total_match / retry + duration_total_build / retry);
}
