#ifndef READ_HAP_H
#define READ_HAP_H
#include <vector>
#include <bits/stdc++.h>
using namespace std;
void parallel_run(int argc, char *argv[]);
void single_run(int argc, char *argv[]);
vector<vector<int> > build_and_match(
    const vector<vector<int> > &haplotypes_sites_matrix,
    int num_threads,
    int L_long_match_min_len);
vector<vector<int> > read_hap(const string &filename);
vector<vector<int> > read_hap_normal(const string &filename);
vector<vector<int> > read_hap_gz(const string &filename);
void print_matches(const vector<vector<int> > &matches);
#endif //READ_HAP_H
