#include <vector>
#include <numeric>
#include <algorithm>

#include "util.h"
using namespace std;

void algorithm_2_BuildPrefixAndDivergenceArrays(const vector<int> &x_k, int k, vector<int> &a, vector<int> &b,
                                                vector<int> &d, vector<int> &e, int M_haps) {
    int u = 0, v = 0;
    int p_val = k + 1, q_val = k + 1;

    for (int i = 0; i < M_haps; ++i) {
        if (d[i] > p_val) p_val = d[i];
        if (d[i] > q_val) q_val = d[i];
        if (x_k[a[i]] == 0) {
            a[u] = a[i];
            d[u] = p_val;
            u++;
            p_val = 0;
        } else {
            b[v] = a[i];
            e[v] = q_val;
            v++;
            q_val = 0;
        }
    }
    copy(b.begin(), b.begin() + v, a.begin() + u);
    copy(e.begin(), e.begin() + v, d.begin() + u);
}

void algorithm_3_ReportLongMatches(
    const vector<int> &x_k_val,
    int N_total_sites,
    int k_current_site,
    int L_min_len,
    const vector<int> &a_arr,
    const vector<int> &d_arr,
    int &i0_val,
    vector<vector<int> > &matches_output_ref,
    int M_total_haps) {
    int u = 0;
    int v = 0;
    int ia = 0;
    int ib = 0;
    int dmin = 0;

    for (int i = 0; i < M_total_haps; ++i) {
        bool condition_met = false;
        int threshold = (k_current_site > L_min_len) ? (k_current_site - L_min_len) : 0;
        if (d_arr[i] > threshold) {
            condition_met = true;
        }

        if (condition_met) {
            if (u > 0 && v > 0) {
                for (ia = i0_val; ia < i; ++ia) {
                    dmin = 0;
                    for (ib = ia + 1; ib < i; ++ib) {
                        if (d_arr[ib] > dmin) dmin = d_arr[ib];
                        if (x_k_val[a_arr[ib]] != x_k_val[a_arr[ia]]) {
                            if (k_current_site > dmin && (k_current_site - dmin >= L_min_len)) {
                                int hap1 = min(a_arr[ia], a_arr[ib]);
                                int hap2 = max(a_arr[ia], a_arr[ib]);
                                matches_output_ref.push_back({k_current_site, hap1, hap2, (k_current_site - dmin)});
                            }
                        }
                    }
                }
            }
            u = 0;
            v = 0;
            i0_val = i;
        }
        if (x_k_val[a_arr[i]] == 0) {
            u++;
        } else {
            v++;
        }
    }

    if (u > 0 && v > 0) {
        for (ia = i0_val; ia < M_total_haps; ++ia) {
            dmin = 0;
            for (ib = ia + 1; ib < M_total_haps; ++ib) {
                if (d_arr[ib] > dmin) dmin = d_arr[ib];
                if (x_k_val[a_arr[ib]] != x_k_val[a_arr[ia]]) {
                    if (k_current_site > dmin && (k_current_site - dmin >= L_min_len)) {
                        int hap1 = min(a_arr[ia], a_arr[ib]);
                        int hap2 = max(a_arr[ia], a_arr[ib]);
                        matches_output_ref.push_back({k_current_site, hap1, hap2, (k_current_site - dmin)});
                    }
                }
            }
        }
    }
}

vector<vector<int> > run_pbwt_corely(
    const vector<vector<int> > &Xt,
    int k_start_site,
    int k_stop_site,
    vector<int> initial_a,
    vector<int> initial_d,
    int M_haps,
    int N_sites,
    int L_min_len,
    bool report_matches_flag,
    vector<vector<int> > &matches_collector) {
    vector<int> a(M_haps), b(M_haps);
    if (!initial_a.empty()) {
        a = initial_a;
    } else {
        a.resize(M_haps);
        iota(a.begin(), a.end(), 0);
    }

    vector<int> d(M_haps), e(M_haps);
    if (!initial_d.empty()) {
        d = initial_d;
    } else {
        d.assign(M_haps, k_start_site);
    }

    int i0 = k_start_site > 0 ? M_haps : 0;

    for (int k = k_start_site; k < k_stop_site; ++k) {
        if (report_matches_flag) {
            algorithm_3_ReportLongMatches(Xt[k], N_sites, k, L_min_len, a, d, i0, matches_collector, M_haps);
        }
        algorithm_2_BuildPrefixAndDivergenceArrays(Xt[k], k, a, b, d, e, M_haps);
    }

    if (report_matches_flag && k_stop_site == N_sites && k_stop_site > 0) {
        algorithm_3_ReportLongMatches(Xt[N_sites - 1], N_sites, N_sites, L_min_len, a, d, i0,
                                                 matches_collector, M_haps);
    }
    return {a, d};
}

vector<vector<int> > build_and_match(const vector<vector<int> > &hap_map_original, int num_threads, int L) {
    if (hap_map_original.empty() || hap_map_original[0].empty() || L <= 0) {
        return {};
    }

    int M = hap_map_original.size();
    int N = hap_map_original[0].size();

    vector<vector<int> > Xt(N, vector<int>(M));
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < N; ++c) {
            Xt[c][r] = hap_map_original[r][c];
        }
    }

    vector<vector<int> > all_final_matches;

    vector<int> a_initial(M);
    iota(a_initial.begin(), a_initial.end(), 0);
    vector<int> d_initial(M, 0);

    run_pbwt_corely(Xt, 0, N, a_initial, d_initial, M, N, L, true, all_final_matches);

    if (!all_final_matches.empty()) {
        sort(all_final_matches.begin(), all_final_matches.end());
        all_final_matches.erase(unique(all_final_matches.begin(), all_final_matches.end()), all_final_matches.end());
    }

    return all_final_matches;
}

int main(int argc, char *argv[]) {
    single_run(argc, argv);
}
