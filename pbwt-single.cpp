#include <vector>
#include <numeric>
#include <algorithm>
#include <memory>
#include <cstring>
#include <cstdio>
#include "util.h"
using namespace std;

void algorithm_2_BuildPrefixAndDivergenceArrays(
    const vector<int>& x_k,
    int k,
    vector<int>& a,
    vector<int>& b,
    vector<int>& d,
    vector<int>& e,
    int M_haps)
{
    int u = 0, v = 0;
    int p_val = k + 1, q_val = k + 1;

    // Unroll small loops for better CPU pipelining
    const int step = 4;
    const int main_loop_limit = M_haps - (M_haps % step);

    // Main loop with unrolling for better instruction-level parallelism
    for (int i = 0; i < main_loop_limit; i += step) {
        // Prefetch data for better cache utilization
        __builtin_prefetch(&a[i + step], 0);
        __builtin_prefetch(&d[i + step], 0);
        __builtin_prefetch(&x_k[a[i + step]], 0);

        // Process 4 elements at once
        for (int j = 0; j < step; ++j) {
            const int idx = i + j;
            const int current_a = a[idx];
            const int current_d = d[idx];

            // Update max values
            p_val = max(p_val, current_d);
            q_val = max(q_val, current_d);

            // Split based on x_k value (0 or 1)
            if (x_k[current_a] == 0) {
                a[u] = current_a;
                d[u] = p_val;
                ++u;
                p_val = 0;
            } else {
                b[v] = current_a;
                e[v] = q_val;
                ++v;
                q_val = 0;
            }
        }
    }

    // Handle remaining elements
    for (int i = main_loop_limit; i < M_haps; ++i) {
        const int current_a = a[i];
        const int current_d = d[i];

        p_val = max(p_val, current_d);
        q_val = max(q_val, current_d);

        if (x_k[current_a] == 0) {
            a[u] = current_a;
            d[u] = p_val;
            ++u;
            p_val = 0;
        } else {
            b[v] = current_a;
            e[v] = q_val;
            ++v;
            q_val = 0;
        }
    }

    // Use memcpy instead of std::copy for performance
    if (v > 0) {
        memcpy(a.data() + u, b.data(), v * sizeof(int));
        memcpy(d.data() + u, e.data(), v * sizeof(int));
    }
}

void algorithm_3_ReportLongMatches(
    const vector<int>& x_k_val,
    int N_total_sites,
    int k_current_site,
    int L_min_len,
    const vector<int>& a_arr,
    const vector<int>& d_arr,
    int& i0_val,
    vector<vector<int>>& matches_output_ref,
    int M_total_haps)
{
    // Pre-calculate threshold once
    const int threshold = max(0, k_current_site - L_min_len);

    // Pre-allocate more capacity for matches to reduce reallocations
    const int estimated_new_matches = (M_total_haps * M_total_haps) / 100;
    if ((int)matches_output_ref.capacity() < (int)matches_output_ref.size() + estimated_new_matches) {
        matches_output_ref.reserve(matches_output_ref.size() + estimated_new_matches);
    }

    int u = 0;
    int v = 0;

    // Main loop with optimized condition checking
    for (int i = 0; i < M_total_haps; ++i) {
        const bool condition_met = (d_arr[i] > threshold);

        if (condition_met) {
            if (u > 0 && v > 0) {
                // Optimize inner loops using early breaks and better cache patterns
                for (int ia = i0_val; ia < i; ++ia) {
                    int dmin = 0;
                    const int hap_ia = a_arr[ia];
                    const int allele_ia = x_k_val[hap_ia];

                    for (int ib = ia + 1; ib < i; ++ib) {
                        dmin = max(dmin, d_arr[ib]);

                        // Skip unnecessary comparisons
                        if (k_current_site <= dmin || (k_current_site - dmin < L_min_len)) {
                            continue;
                        }

                        const int hap_ib = a_arr[ib];
                        if (x_k_val[hap_ib] != allele_ia) {
                            matches_output_ref.push_back({
                                k_current_site,
                                min(hap_ia, hap_ib),
                                max(hap_ia, hap_ib),
                                k_current_site - dmin
                            });
                        }
                    }
                }
            }
            u = 0;
            v = 0;
            i0_val = i;
        }

        // Update counters
        if (x_k_val[a_arr[i]] == 0) {
            ++u;
        } else {
            ++v;
        }
    }

    // Process final block
    if (u > 0 && v > 0) {
        for (int ia = i0_val; ia < M_total_haps; ++ia) {
            int dmin = 0;
            const int hap_ia = a_arr[ia];
            const int allele_ia = x_k_val[hap_ia];

            for (int ib = ia + 1; ib < M_total_haps; ++ib) {
                dmin = max(dmin, d_arr[ib]);

                // Skip unnecessary comparisons
                if (k_current_site <= dmin || (k_current_site - dmin < L_min_len)) {
                    continue;
                }

                const int hap_ib = a_arr[ib];
                if (x_k_val[hap_ib] != allele_ia) {
                    matches_output_ref.push_back({
                        k_current_site,
                        min(hap_ia, hap_ib),
                        max(hap_ia, hap_ib),
                        k_current_site - dmin
                    });
                }
            }
        }
    }
}

vector<vector<int>> run_pbwt_corely(
    const vector<vector<int>>& Xt,
    int k_start_site,
    int k_stop_site,
    vector<int> initial_a,
    vector<int> initial_d,
    int M_haps,
    int N_sites,
    int L_min_len,
    bool report_matches_flag,
    vector<vector<int>>& matches_collector)
{
    // Preallocate vectors with exact sizes
    vector<int> a(M_haps), b(M_haps), d(M_haps), e(M_haps);

    // Initialize arrays
    if (!initial_a.empty()) {
        a = move(initial_a);  // Use move instead of copy
    } else {
        iota(a.begin(), a.end(), 0);
    }

    if (!initial_d.empty()) {
        d = move(initial_d);  // Use move instead of copy
    } else {
        fill(d.begin(), d.end(), k_start_site);  // fill is faster than assign
    }

    int i0 = k_start_site > 0 ? M_haps : 0;

    // Reserve memory for matches to reduce reallocations
    if (report_matches_flag) {
        const int estimated_matches = M_haps * (k_stop_site - k_start_site) / 20;
        matches_collector.reserve(matches_collector.size() + estimated_matches);
    }

    // Main PBWT loop
    for (int k = k_start_site; k < k_stop_site; ++k) {
        if (report_matches_flag) {
            algorithm_3_ReportLongMatches(Xt[k], N_sites, k, L_min_len, a, d, i0, matches_collector, M_haps);
        }
        algorithm_2_BuildPrefixAndDivergenceArrays(Xt[k], k, a, b, d, e, M_haps);
    }

    // Handle final report
    if (report_matches_flag && k_stop_site == N_sites && k_stop_site > 0) {
        algorithm_3_ReportLongMatches(Xt[N_sites - 1], N_sites, N_sites, L_min_len, a, d, i0,
                                      matches_collector, M_haps);
    }

    return {a, d};
}

vector<vector<int>> build_and_match(const vector<vector<int>>& hap_map_original, int num_threads, int L) {
    // Early validation
    if (hap_map_original.empty() || hap_map_original[0].empty() || L <= 0) {
        return {};
    }

    const int M = hap_map_original.size();
    const int N = hap_map_original[0].size();

    // Optimize the transpose operation for better cache utilization
    vector<vector<int>> Xt(N, vector<int>(M));
    for (int c = 0; c < N; ++c) {
        auto& row = Xt[c];
        for (int r = 0; r < M; ++r) {
            row[r] = hap_map_original[r][c];
        }
    }

    // Pre-allocate vectors with exact sizes
    vector<int> a_initial(M);
    iota(a_initial.begin(), a_initial.end(), 0);

    vector<int> d_initial(M, 0);

    // Optimize memory allocation for matches
    vector<vector<int>> all_final_matches;
    const int estimated_initial_matches = min(M * N / 20, 10000);
    all_final_matches.reserve(estimated_initial_matches);

    // Run the PBWT algorithm
    run_pbwt_corely(Xt, 0, N, a_initial, d_initial, M, N, L, true, all_final_matches);

    // Optimize post-processing
    if (!all_final_matches.empty()) {
        // Sort and remove duplicates
        sort(all_final_matches.begin(), all_final_matches.end());
        all_final_matches.erase(unique(all_final_matches.begin(), all_final_matches.end()), all_final_matches.end());
    }

    return all_final_matches;
}

std::vector<std::vector<int>> algorithm_4_ReportMaximalMatches(
    const std::vector<std::vector<int>>& hap_matrix,
    const std::vector<int>& p_k,
    const std::vector<int>& d_k,
    int k
) {
    int M = hap_matrix.size();
    int N = hap_matrix[0].size();
    std::vector<std::vector<int>> matches;

    // sentinel at boundaries
    std::vector<int> d_ext(M + 2, k + 1);
    for (int i = 0; i < M; ++i) d_ext[i + 1] = d_k[i];

    for (int i = 0; i < M; ++i) {
        int di = d_ext[i + 1];
        int dip1 = d_ext[i + 2];
        int yi = hap_matrix[p_k[i]][k];
        int m = i - 1, n = i + 1;

        // ---- scan down ----
        bool skip = false;
        if (di <= dip1) {
            while (m >= 0 && d_ext[m + 1] <= di) {
                int ym = hap_matrix[p_k[m]][k];
                if (ym == yi && k != N - 1) {
                    skip = true;
                    break;
                }
                m--;
            }
        }
        // ---- scan up ----
        if (!skip && di >= dip1) {
            while (n < M && d_ext[n + 1] <= dip1) {
                int yn = hap_matrix[p_k[n]][k];
                if (yn == yi && k != N - 1) {
                    skip = true;
                    break;
                }
                n++;
            }
        }
        if (skip) continue;

        // report matches
        for (int j = m + 1; j < i; ++j) {
            int h1 = std::min(p_k[i], p_k[j]);
            int h2 = std::max(p_k[i], p_k[j]);
            int length = k - di;
            if (length > 0) matches.push_back({di, h1, h2, length});
        }
        for (int j = i + 1; j < n; ++j) {
            int h1 = std::min(p_k[i], p_k[j]);
            int h2 = std::max(p_k[i], p_k[j]);
            int length = k - dip1;
            if (length > 0) matches.push_back({dip1, h1, h2, length});
        }
    }
    return matches;
}




std::vector<std::vector<int>> build_and_match_maximal(
    const std::vector<std::vector<int>>& hap_map_original,
    int num_threads,
    int L
) {
    if (hap_map_original.empty() || hap_map_original[0].empty()) {
        std::cerr << "ERROR: Input haplotype matrix is empty.\n";
        return {};
    }

    const int M = hap_map_original.size();
    const int N = hap_map_original[0].size();

    for (int i = 1; i < M; ++i) {
        if ((int)hap_map_original[i].size() != N) {
            std::cerr << "ERROR: Inconsistent row length in haplotype matrix at row " << i << ".\n";
            return {};
        }
    }

    std::vector<std::vector<int>> Xt(N, std::vector<int>(M));
    for (int c = 0; c < N; ++c)
        for (int r = 0; r < M; ++r)
            Xt[c][r] = hap_map_original[r][c];

    std::vector<int> a(M); std::iota(a.begin(), a.end(), 0);
    std::vector<int> b(M), d(M, 0), e(M);

    std::vector<std::vector<int>> maximal_matches;
    const int estimated_initial_matches = std::min(M * N / 10, 10000);
    maximal_matches.reserve(estimated_initial_matches);

    for (int k = 0; k < N; ++k) {
        auto k_matches = algorithm_4_ReportMaximalMatches(hap_map_original, a, d, k);
        maximal_matches.insert(maximal_matches.end(), k_matches.begin(), k_matches.end());

        algorithm_2_BuildPrefixAndDivergenceArrays(Xt[k], k, a, b, d, e, M);
    }

    if (!maximal_matches.empty()) {
        std::sort(maximal_matches.begin(), maximal_matches.end());
        maximal_matches.erase(std::unique(maximal_matches.begin(), maximal_matches.end()), maximal_matches.end());
    }

    return maximal_matches;
}






int main(int argc, char *argv[]) {
    single_run(argc, argv);
}