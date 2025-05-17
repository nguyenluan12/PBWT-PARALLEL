#ifndef READ_HAP_H
#define READ_HAP_H
#include <vector>
#include <bits/stdc++.h>
using namespace std;
vector<vector<int> > read_hap(const string &filename);
vector<vector<int> > read_hap_normal(const string &filename);
vector<vector<int> > read_hap_gz(const string &filename);
void print_matches(const vector<vector<int> > &matches);
#endif //READ_HAP_H
