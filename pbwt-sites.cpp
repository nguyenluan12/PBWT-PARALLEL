#include <vector>
#include <numeric>
#include <algorithm>
#include <omp.h>
#include <string>
#include <utility> // For std::pair

#include "util.h"

using namespace std;
using namespace chrono;

void algorithm_2_BuildPrefixAndDivergenceArrays(
    const vector<int> &x_k,
    int k, vector<int> &a,
    vector<int> &b,
    vector<int> &d,
    vector<int> &e,
    int M_haps) {
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
#pragma omp critical
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
#pragma omp critical
                        matches_output_ref.push_back({k_current_site, hap1, hap2, (k_current_site - dmin)});
                    }
                }
            }
        }
    }
}

void fill_rppa(vector<int> &rppa, const vector<int> &ppa) {
    const int N_rppa = ppa.size();
    for (int i = 0; i < N_rppa; ++i) {
        rppa[ppa[i]] = i;
    }
}

pair<vector<int>, vector<int> > run_pbwt_sequentially(
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
        algorithm_3_ReportLongMatches(Xt[N_sites - 1], N_sites, N_sites, L_min_len, a, d, i0, matches_collector,
                                      M_haps);
    }
    return {a, d};
}

void fix_a_d_range_internal(
    int group_start_idx, int group_stop_idx,
    const std::vector<int> &prev_a_state,
    const std::vector<int> &prev_d_state,
    const std::vector<int> &rppa_lookup,
    std::vector<int> &current_a_to_fix,
    std::vector<int> &current_d_to_fix)
{
    int sz = group_stop_idx - group_start_idx;
    std::vector<int> prev_pos(sz);
    for (int j = 0; j < sz; ++j) {
        prev_pos[j] = rppa_lookup[current_a_to_fix[group_start_idx + j]];
    }
    std::sort(prev_pos.begin(), prev_pos.end());

    // Sắp lại thứ tự của a theo thứ tự prev_a_state
    for (int j = 0; j < sz; ++j) {
        current_a_to_fix[group_start_idx + j] = prev_a_state[prev_pos[j]];
    }

    // Cập nhật lại d theo max của prev_d_state
    for (int j = 1; j < sz; ++j) {
        int scan_start = prev_pos[j - 1] + 1;
        int scan_stop  = prev_pos[j] + 1;
        if (scan_start < scan_stop && scan_stop <= (int)prev_d_state.size()) {
            current_d_to_fix[group_start_idx + j] = *std::max_element(
                prev_d_state.begin() + scan_start, prev_d_state.begin() + scan_stop
            );
        } else if (scan_start == scan_stop && scan_start > 0 && scan_start <= (int)prev_d_state.size()) {
            current_d_to_fix[group_start_idx + j] = prev_d_state[scan_start - 1];
        }
    }
}


void fix_checkpoints(
    std::vector<std::vector<int>> &checkpoint_as,
    std::vector<std::vector<int>> &checkpoint_ds,
    const std::vector<int> &checkpoint_k_values,
    int M_haps)
{
    if (checkpoint_as.empty()) return;
    std::vector<int> rppa(M_haps);

    // Sửa từ checkpoint thứ 2 trở đi
    for (int i = 1; i < (int)checkpoint_as.size(); ++i) {
        std::vector<int> &a_to_fix = checkpoint_as[i];
        std::vector<int> &d_to_fix = checkpoint_ds[i];
        const std::vector<int> &prev_a = checkpoint_as[i - 1];
        const std::vector<int> &prev_d = checkpoint_ds[i - 1];
        int k_val_of_prev_state = checkpoint_k_values[i - 1];

        // Xây dựng rppa (inverse of prev_a)
        for (int idx = 0; idx < M_haps; ++idx)
            rppa[prev_a[idx]] = idx;

        int group_start = 0;
        for (int j = 0; j < M_haps; ++j) {
            if (d_to_fix[j] != k_val_of_prev_state) {
                int group_size = j - group_start;
                if (group_size > 1) {
                    fix_a_d_range_internal(group_start, j, prev_a, prev_d, rppa, a_to_fix, d_to_fix);
                }
                group_start = j;
            }
        }
        // Sửa nhóm cuối cùng nếu cần
        if (M_haps - group_start > 1) {
            fix_a_d_range_internal(group_start, M_haps, prev_a, prev_d, rppa, a_to_fix, d_to_fix);
        }
    }
}


vector<vector<int> > build_and_match(const vector<vector<int> > &hap_map_original, int num_threads, int L) {
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);

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
    vector<int> checkpoint_positions_vec;

    if (num_threads > 1 && N > 0) {
        int step = N / num_threads;
        if (step == 0 && N > 0) step = N; // Ensure at least one segment or full if N < num_threads
        for (int i = 1; i < num_threads; ++i) {
            int pos = step * i;
            if (pos < N && pos > 0) {
                checkpoint_positions_vec.push_back(pos);
            } else {
                break;
            }
        }
        sort(checkpoint_positions_vec.begin(), checkpoint_positions_vec.end());
        checkpoint_positions_vec.erase(unique(checkpoint_positions_vec.begin(), checkpoint_positions_vec.end()),
                                       checkpoint_positions_vec.end());
    }

    vector<vector<int> > checkpoint_as_vec(checkpoint_positions_vec.size());
    vector<vector<int> > checkpoint_ds_vec(checkpoint_positions_vec.size());

    if (num_threads > 1 && !checkpoint_positions_vec.empty()) {
#pragma omp parallel for
        for (int i = 0; i < (int) checkpoint_positions_vec.size(); ++i) {
            int k_start = (i == 0) ? 0 : checkpoint_positions_vec[i - 1];
            int k_stop = checkpoint_positions_vec[i];

            vector<int> seg_a_init(M);
            iota(seg_a_init.begin(), seg_a_init.end(), 0);
            vector<int> seg_d_init(M, k_start);

            if (k_start < k_stop) {
                pair<vector<int>, vector<int> > ad_pair = run_pbwt_sequentially(
                    Xt, k_start, k_stop, seg_a_init, seg_d_init, M, N, L, false, all_final_matches);
                checkpoint_as_vec[i] = ad_pair.first;
                checkpoint_ds_vec[i] = ad_pair.second;
            }
        }
        fix_checkpoints(checkpoint_as_vec, checkpoint_ds_vec, checkpoint_positions_vec, M);
    }

    vector<vector<vector<int> > > parallel_match_results(num_threads);
#pragma omp parallel for
    for (int i = 0; i < num_threads; ++i) {
        int k_start_segment, k_stop_segment;
        vector<int> a_segment_initial;
        vector<int> d_segment_initial;

        if (num_threads == 1 || checkpoint_positions_vec.empty()) {
            k_start_segment = 0;
            k_stop_segment = N;
            a_segment_initial.resize(M);
            iota(a_segment_initial.begin(), a_segment_initial.end(), 0);
            d_segment_initial.assign(M, 0);
        } else {
            k_start_segment = (i == 0) ? 0 : checkpoint_positions_vec[i - 1];
            k_stop_segment = (i == (int) checkpoint_positions_vec.size()) ? N : checkpoint_positions_vec[i];

            if (i == 0) {
                a_segment_initial.resize(M);
                iota(a_segment_initial.begin(), a_segment_initial.end(), 0);
                d_segment_initial.assign(M, 0);
            } else {
                if (i - 1 < (int) checkpoint_as_vec.size()) {
                    // Ensure checkpoint exists
                    a_segment_initial = checkpoint_as_vec[i - 1];
                    d_segment_initial = checkpoint_ds_vec[i - 1];
                } else {
                    // Fallback if checkpoint logic had issues, process from start (less efficient but safer)
                    a_segment_initial.resize(M);
                    iota(a_segment_initial.begin(), a_segment_initial.end(), 0);
                    d_segment_initial.assign(M, k_start_segment); // d relative to current start
                }
            }
        }
        if (k_start_segment < k_stop_segment) {
            run_pbwt_sequentially(Xt, k_start_segment, k_stop_segment, a_segment_initial, d_segment_initial, M, N, L,
                                  true, parallel_match_results[i]);
        }
    }

    for (const auto &thread_matches: parallel_match_results) {
        all_final_matches.insert(all_final_matches.end(), thread_matches.begin(), thread_matches.end());
    }

    if (!all_final_matches.empty()) {
        sort(all_final_matches.begin(), all_final_matches.end());
        all_final_matches.erase(unique(all_final_matches.begin(), all_final_matches.end()), all_final_matches.end());
    }

    return all_final_matches;
}
// Parallel PBWT - support for both long and set-maximal matches


void algorithm_4_ReportSetMaximalMatches(
    const std::vector<int> &x_k_val,
    int N_total_sites,
    int k_current_site,
    int L_min_len,
    const std::vector<int> &a_arr,
    const std::vector<int> &d_arr,
    std::vector<std::vector<int>> &matches_output_ref,
    int M_total_haps) 
{
    // Bổ sung sentinel ở hai đầu mảng d để thuận tiện cho quét m/n.
    std::vector<int> d(M_total_haps + 2, k_current_site + 1);
    for (int i = 0; i < M_total_haps; ++i)
        d[i + 1] = d_arr[i];

    for (int i = 0; i < M_total_haps; ++i) {
        int m = i - 1;
        int n = i + 1;
        bool skip_match = false;

        // Quét xuống (scan down the array)
        if (d[i + 1] <= d[i + 2]) {
            while (m >= 0 && d[m + 1] <= d[i + 1]) {
                // Nếu gặp genotype giống và chưa phải cuối, bỏ qua vị trí i
                if (x_k_val[a_arr[m]] == x_k_val[a_arr[i]] && k_current_site != N_total_sites - 1) {
                    skip_match = true;
                    break;
                }
                --m;
            }
        }

        // Quét lên (scan up the array)
        if (!skip_match && d[i + 1] >= d[i + 2]) {
            while (n < M_total_haps && d[n + 1] <= d[i + 2]) {
                if (x_k_val[a_arr[n]] == x_k_val[a_arr[i]] && k_current_site != N_total_sites - 1) {
                    skip_match = true;
                    break;
                }
                ++n;
            }
        }

        if (skip_match) continue;

        // Báo cáo các match từ m+1 đến i-1 với i (theo chiều xuống)
        for (int j = m + 1; j < i; ++j) {
            int match_start = d[i + 1];
            int match_len = k_current_site - match_start;
            if (match_len >= L_min_len) {
                int h1 = std::min(a_arr[i], a_arr[j]);
                int h2 = std::max(a_arr[i], a_arr[j]);
#pragma omp critical
                matches_output_ref.push_back({match_start, h1, h2, match_len});
            }
        }

        // Báo cáo các match từ i+1 đến n-1 với i (theo chiều lên)
        for (int j = i + 1; j < n; ++j) {
            int match_start = d[i + 2];
            int match_len = k_current_site - match_start;
            if (match_len >= L_min_len) {
                int h1 = std::min(a_arr[i], a_arr[j]);
                int h2 = std::max(a_arr[i], a_arr[j]);
#pragma omp critical
                matches_output_ref.push_back({match_start, h1, h2, match_len});
            }
        }
    }
}


// Sequential PBWT run specifically for set-maximal matches
pair<vector<int>, vector<int>> run_pbwt_sequentially_setmax(
    const vector<vector<int>> &Xt,
    int k_start_site,
    int k_stop_site,
    vector<int> initial_a,
    vector<int> initial_d,
    int M_haps,
    int N_sites,
    int L_min_len,
    bool report_matches_flag,
    vector<vector<int>> &matches_collector) {
    
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

    for (int k = k_start_site; k < k_stop_site; ++k) {
        if (report_matches_flag) {
            algorithm_4_ReportSetMaximalMatches(Xt[k], N_sites, k, L_min_len, a, d, matches_collector, M_haps);
        }
        algorithm_2_BuildPrefixAndDivergenceArrays(Xt[k], k, a, b, d, e, M_haps);
    }

    // Handle final position if needed
    if (report_matches_flag && k_stop_site == N_sites && k_stop_site > 0) {
        algorithm_4_ReportSetMaximalMatches(Xt[N_sites - 1], N_sites, N_sites, L_min_len, a, d, matches_collector, M_haps);
    }
    
    return {a, d};
}

std::vector<std::vector<int>> build_and_match_maximal(const std::vector<std::vector<int>>& hap_map_original,
                                                         int num_threads,
                                                         int L_min_len) {
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);

    if (hap_map_original.empty() || hap_map_original[0].empty() || L_min_len <= 0) {
        return {};
    }

    int M = hap_map_original.size();
    int N = hap_map_original[0].size();

    // Transpose matrix to have sites as rows and haplotypes as columns
    std::vector<std::vector<int>> Xt(N, std::vector<int>(M));
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < N; ++c) {
            Xt[c][r] = hap_map_original[r][c];
        }
    }

    std::vector<std::vector<int>> all_final_matches;
    std::vector<int> checkpoint_positions_vec;

    // Setup checkpoints for parallel processing
    if (num_threads > 1 && N > 0) {
        int step = N / num_threads;
        if (step == 0 && N > 0) step = N;
        for (int i = 1; i < num_threads; ++i) {
            int pos = step * i;
            if (pos < N && pos > 0) {
                checkpoint_positions_vec.push_back(pos);
            } else {
                break;
            }
        }
        std::sort(checkpoint_positions_vec.begin(), checkpoint_positions_vec.end());
        checkpoint_positions_vec.erase(std::unique(checkpoint_positions_vec.begin(), checkpoint_positions_vec.end()),
                                       checkpoint_positions_vec.end());
    }

    std::vector<std::vector<int>> checkpoint_as_vec(checkpoint_positions_vec.size());
    std::vector<std::vector<int>> checkpoint_ds_vec(checkpoint_positions_vec.size());

    // Phase A: Generate approximate a and d arrays in parallel
    if (num_threads > 1 && !checkpoint_positions_vec.empty()) {
#pragma omp parallel for
        for (int i = 0; i < (int) checkpoint_positions_vec.size(); ++i) {
            int k_start = (i == 0) ? 0 : checkpoint_positions_vec[i - 1];
            int k_stop = checkpoint_positions_vec[i];

            std::vector<int> seg_a_init(M);
            std::iota(seg_a_init.begin(), seg_a_init.end(), 0);
            std::vector<int> seg_d_init(M, k_start);

            if (k_start < k_stop) {
                std::vector<std::vector<int>> dummy_matches;
                auto ad_pair = run_pbwt_sequentially_setmax(
                    Xt, k_start, k_stop, seg_a_init, seg_d_init,
                    M, N, L_min_len, false, dummy_matches);
                checkpoint_as_vec[i] = ad_pair.first;
                checkpoint_ds_vec[i] = ad_pair.second;
            }
        }
        
        // Phase B: Sequential correction of approximate arrays
        fix_checkpoints(checkpoint_as_vec, checkpoint_ds_vec, checkpoint_positions_vec, M);
    }

    // Phase C: Parallel execution with corrected arrays
    std::vector<std::vector<std::vector<int>>> parallel_match_results(num_threads);

#pragma omp parallel for
    for (int i = 0; i < num_threads; ++i) {
        int k_start_segment, k_stop_segment;
        std::vector<int> a_segment_initial, d_segment_initial;

        if (num_threads == 1 || checkpoint_positions_vec.empty()) {
            // Single thread case
            k_start_segment = 0;
            k_stop_segment = N;
            a_segment_initial.resize(M);
            std::iota(a_segment_initial.begin(), a_segment_initial.end(), 0);
            d_segment_initial.assign(M, 0);
        } else {
            // Multi-thread case
            k_start_segment = (i == 0) ? 0 : checkpoint_positions_vec[i - 1];
            k_stop_segment = (i == (int) checkpoint_positions_vec.size()) ? N : checkpoint_positions_vec[i];

            if (i == 0) {
                a_segment_initial.resize(M);
                std::iota(a_segment_initial.begin(), a_segment_initial.end(), 0);
                d_segment_initial.assign(M, 0);
            } else {
                if (i - 1 < (int) checkpoint_as_vec.size()) {
                    a_segment_initial = checkpoint_as_vec[i - 1];
                    d_segment_initial = checkpoint_ds_vec[i - 1];
                } else {
                    // Fallback
                    a_segment_initial.resize(M);
                    std::iota(a_segment_initial.begin(), a_segment_initial.end(), 0);
                    d_segment_initial.assign(M, k_start_segment);
                }
            }
        }

        if (k_start_segment < k_stop_segment) {
            run_pbwt_sequentially_setmax(Xt, k_start_segment, k_stop_segment, 
                                       a_segment_initial, d_segment_initial, 
                                       M, N, L_min_len, true, parallel_match_results[i]);
        }
    }

    // Combine results from all threads
    for (const auto& thread_matches : parallel_match_results) {
        all_final_matches.insert(all_final_matches.end(), thread_matches.begin(), thread_matches.end());
    }

    // Remove duplicates and sort
    if (!all_final_matches.empty()) {
        std::sort(all_final_matches.begin(), all_final_matches.end());
        all_final_matches.erase(std::unique(all_final_matches.begin(), all_final_matches.end()), all_final_matches.end());
    }

    return all_final_matches;
}


int main(int argc, char *argv[]) {
    parallel_run(argc, argv);
}
