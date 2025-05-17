#ifndef PBWT_HAPS_H
#define PBWT_HAPS_H
#include <vector>

using namespace std;

pair<vector<int>, vector<int>> parallel_prefix_array_computation(
    const vector<int>& ppa_at_prev_site,
    const vector<vector<int>>& data_matrix,
    int site_idx,
    int N,
    int num_threads);

vector<int> parallel_divergence_array_computation(
    const vector<int>& ppa_at_prev_site,
    const vector<int>& mlens_at_prev_site,
    const vector<int>& ps_holder_for_ppa_prev_site,
    const vector<vector<int>>& data_matrix,
    int site_idx,
    int N,
    int num_threads);
#endif //PBWT_HAPS_H
