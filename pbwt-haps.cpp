#include <vector>
#include <numeric>
#include <algorithm>
#include <omp.h>
#include <cmath>
#include <limits>
#include <cstring> // For memcpy

#include "util.h"

using namespace std;

// Forward declarations
void algorithm_2_BuildPrefixAndDivergenceArrays(
    const std::vector<int>& x_k,
    int k,
    std::vector<int>& a,
    std::vector<int>& b,
    std::vector<int>& d,
    std::vector<int>& e,
    int M_haps);

std::vector<std::vector<int>> algorithm_4_ReportMaximalMatches(
    const std::vector<std::vector<int>>& hap_matrix,
    const std::vector<int>& p_k,
    const std::vector<int>& d_k,
    int k
);

// Implementation of Algorithm 2 from the paper
void algorithm_2_BuildPrefixAndDivergenceArrays(
    const std::vector<int>& x_k,
    int k,
    std::vector<int>& a,
    std::vector<int>& b,
    std::vector<int>& d,
    std::vector<int>& e,
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
            p_val = std::max(p_val, current_d);
            q_val = std::max(q_val, current_d);

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

        p_val = std::max(p_val, current_d);
        q_val = std::max(q_val, current_d);

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


void core_p_arr(
    const vector<vector<int>>& haplotypes_sites_matrix,
    int site_idx,
    int n_hap,
    int num_threads,
    int block_size,
    const vector<int>& old_pda_parr,
    vector<int>& new_pda_parr,
    int& new_pda_zerocnt,
    vector<int>& ps_holder,
    vector<int>& offsets_holder)
{
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 0; i < num_threads; ++i) {
        int s = i * block_size;
        if (s >= n_hap) continue;
        int e = min(s + block_size, n_hap);

        if (s < e) {
            const int hap_s = old_pda_parr[s];
            ps_holder[s] = (haplotypes_sites_matrix[hap_s][site_idx] == 0) ? 0 : 1;

            for (int k = s + 1; k < e; ++k) {
                const int hap_k = old_pda_parr[k];
                ps_holder[k] = ps_holder[k - 1] + (haplotypes_sites_matrix[hap_k][site_idx] != 0);
            }
            offsets_holder[i] = ps_holder[e - 1];
        } else {
            offsets_holder[i] = 0;
        }
    }

    // Prefix sum on offsets
    for (int i = 1; i < num_threads; ++i) {
        offsets_holder[i] += offsets_holder[i - 1];
    }

    // Adjust ps values for each block
    #pragma omp parallel for
    for (int i = 1; i < num_threads; ++i) {
        int s = i * block_size;
        if (s >= n_hap) continue;
        int e = min(s + block_size, n_hap);
        for (int k = s; k < e; ++k) {
            ps_holder[k] += offsets_holder[i - 1];
        }
    }

    // Calculate zerocnt and prepare for reordering
    if (n_hap > 0) {
        new_pda_zerocnt = n_hap - ps_holder[n_hap - 1];
    } else {
        new_pda_zerocnt = 0;
    }
    int one_off = new_pda_zerocnt - 1;

    // Reorder array based on bit values
    #pragma omp parallel for
    for (int i = 0; i < n_hap; ++i) {
        const int hap_i = old_pda_parr[i];
        if (haplotypes_sites_matrix[hap_i][site_idx] == 0) {
            new_pda_parr[i - ps_holder[i]] = hap_i;
        } else {
            new_pda_parr[ps_holder[i] + one_off] = hap_i;
        }
    }
}

void initial_sort(
    const vector<vector<int>>& haplotypes_sites_matrix,
    int site_idx,
    int n_hap,
    int num_threads,
    int block_size,
    vector<int>& old_pda_parr,
    vector<int>& old_pda_mlens,
    int& old_pda_zerocnt,
    vector<int>& ps_holder,
    vector<int>& offsets_holder)
{
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 0; i < num_threads; ++i) {
        int s = i * block_size;
        if (s >= n_hap) continue;
        int e = min(s + block_size, n_hap);

        if (s < e) {
            ps_holder[s] = (haplotypes_sites_matrix[s][site_idx] == 0) ? 0 : 1;

            for (int k = s + 1; k < e; ++k) {
                ps_holder[k] = ps_holder[k - 1] + (haplotypes_sites_matrix[k][site_idx] != 0);
            }
            offsets_holder[i] = ps_holder[e - 1];
        } else {
            offsets_holder[i] = 0;
        }
    }

    // Prefix sum on offsets
    for (int i = 1; i < num_threads; ++i) {
        offsets_holder[i] += offsets_holder[i - 1];
    }

    // Adjust ps values for each block
    #pragma omp parallel for
    for (int i = 1; i < num_threads; ++i) {
        int s = i * block_size;
        if (s >= n_hap) continue;
        int e = min(s + block_size, n_hap);
        for (int k = s; k < e; ++k) {
            ps_holder[k] += offsets_holder[i - 1];
        }
    }

    // Calculate zerocnt and prepare for reordering
    if (n_hap > 0) {
        old_pda_zerocnt = n_hap - ps_holder[n_hap - 1];
    } else {
        old_pda_zerocnt = 0;
    }
    int one_off = old_pda_zerocnt - 1;

    // Initialize arrays
    #pragma omp parallel for
    for (int i = 0; i < n_hap; ++i) {
        if (haplotypes_sites_matrix[i][site_idx] == 0) {
            old_pda_parr[i - ps_holder[i]] = i;
            old_pda_mlens[i] = 1;
        } else {
            old_pda_parr[ps_holder[i] + one_off] = i;
            old_pda_mlens[i] = 1;
        }
    }

    // Initialize boundary values
    if (n_hap > 0) {
        if (!old_pda_parr.empty()) old_pda_mlens[old_pda_parr[0]] = 0;
        if (old_pda_zerocnt != n_hap) {
            if (old_pda_zerocnt >= 0 && old_pda_zerocnt < static_cast<int>(old_pda_parr.size())) {
                old_pda_mlens[old_pda_parr[old_pda_zerocnt]] = 0;
            }
        }
    }
}

void core_d_arr(
    const vector<vector<int>>& haplotypes_sites_matrix,
    int site_idx,
    int n_hap,
    int num_threads,
    int block_size,
    const vector<int>& old_pda_parr,
    const vector<int>& old_pda_mlens,
    vector<int>& new_pda_mlens,
    const vector<int>& ps_holder
)
{
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 0; i < num_threads; ++i) {
        int s = i * block_size;
        if (s >= n_hap) continue;
        int e = min(s + block_size, n_hap);

        int prv_low_m_zero = -1;
        int prv_low_m_one = -1;
        int h_id;

        if (s == 0) {
            for (int k = s; k < e; ++k) {
                h_id = old_pda_parr[k];
                if (haplotypes_sites_matrix[h_id][site_idx] == 0) {
                    new_pda_mlens[h_id] = min(old_pda_mlens[h_id], prv_low_m_one) + 1;
                    prv_low_m_one = numeric_limits<int>::max();
                    prv_low_m_zero = min(prv_low_m_zero, old_pda_mlens[h_id]);
                } else {
                    new_pda_mlens[h_id] = min(old_pda_mlens[h_id], prv_low_m_zero) + 1;
                    prv_low_m_zero = numeric_limits<int>::max();
                    prv_low_m_one = min(prv_low_m_one, old_pda_mlens[h_id]);
                }
            }
        } else {
            prv_low_m_zero = -1;
            prv_low_m_one = -1;

            bool min_zero_search = true;
            bool min_one_search = true;

            if (ps_holder[s - 1] == s){
                prv_low_m_zero = -1;
                min_zero_search = false;
            }
            if (ps_holder[s - 1] == 0) {
                prv_low_m_one = -1;
                min_one_search = false;
            }

            if (min_zero_search) {
                int seek_index = s - 1;
                while (seek_index >= 0 && haplotypes_sites_matrix[old_pda_parr[seek_index]][site_idx] != 0) {
                    seek_index--;
                }
                if (seek_index >= 0) {
                    prv_low_m_zero = old_pda_mlens[old_pda_parr[seek_index]];
                    while (seek_index >= 0 && haplotypes_sites_matrix[old_pda_parr[seek_index]][site_idx] == 0) {
                        prv_low_m_zero = min(old_pda_mlens[old_pda_parr[seek_index]], prv_low_m_zero);
                        seek_index--;
                    }
                } else {
                    prv_low_m_zero = -1;
                }
            }

            if (min_one_search) {
                int seek_index = s - 1;
                while (seek_index >= 0 && haplotypes_sites_matrix[old_pda_parr[seek_index]][site_idx] == 0) {
                    seek_index--;
                }
                if (seek_index >=0) {
                    prv_low_m_one = old_pda_mlens[old_pda_parr[seek_index]];
                    while (seek_index >= 0 && haplotypes_sites_matrix[old_pda_parr[seek_index]][site_idx] != 0) {
                        prv_low_m_one = min(old_pda_mlens[old_pda_parr[seek_index]], prv_low_m_one);
                        seek_index--;
                    }
                } else {
                    prv_low_m_one = -1;
                }
            }

            for (int k = s; k < e; ++k) {
                h_id = old_pda_parr[k];
                if (haplotypes_sites_matrix[h_id][site_idx] == 0) {
                    new_pda_mlens[h_id] = min(old_pda_mlens[h_id], prv_low_m_one) + 1;
                    prv_low_m_one = numeric_limits<int>::max();
                    prv_low_m_zero = min(prv_low_m_zero, old_pda_mlens[h_id]);
                } else {
                    new_pda_mlens[h_id] = min(old_pda_mlens[h_id], prv_low_m_zero) + 1;
                    prv_low_m_zero = numeric_limits<int>::max();
                    prv_low_m_one = min(prv_low_m_one, old_pda_mlens[h_id]);
                }
            }
        }
    }
}

vector<vector<int>> report_long_matches_sl(
    int site_index_0_based,
    const vector<vector<int>>& haplotypes_sites_matrix,
    int site_idx,
    int n_hap,
    int num_threads,
    int llm_len,
    const vector<int>& old_pda_parr,
    const vector<int>& old_pda_mlens)
{
    // Tối ưu bằng cách dự đoán kích thước
    const int max_expected_matches = min(n_hap * n_hap / 4, 100000);
    vector<vector<vector<int>>> per_thread_matches_storage(num_threads);
    for (auto& thread_storage : per_thread_matches_storage) {
        thread_storage.reserve(max_expected_matches / num_threads);
    }

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(dynamic, 32)
        for (int h = 0; h < n_hap - 1; ++h) {
            int min_val = numeric_limits<int>::max();
            const int hap_h = old_pda_parr[h];
            const int allele_h = haplotypes_sites_matrix[hap_h][site_idx];

            for (int i = h + 1; i < n_hap; ++i) {
                const int hap_i = old_pda_parr[i];
                if (old_pda_mlens[hap_i] < llm_len) {
                    break;
                }
                min_val = min(min_val, old_pda_mlens[hap_i]);

                if (haplotypes_sites_matrix[hap_i][site_idx] != allele_h) {
                    per_thread_matches_storage[thread_id].push_back({
                        site_index_0_based + 1,
                        min(hap_h, hap_i),
                        max(hap_h, hap_i),
                        min_val
                    });
                }
            }
        }
    }

    // Tính toán tổng số matches để dự đoán kích thước
    size_t total_matches = 0;
    for (const auto& thread_list : per_thread_matches_storage) {
        total_matches += thread_list.size();
    }

    vector<vector<int>> aggregated_function_matches;
    aggregated_function_matches.reserve(total_matches);

    // Sử dụng move iterator để tối ưu việc nối
    for (auto& thread_list : per_thread_matches_storage) {
        if (!thread_list.empty()) {
            aggregated_function_matches.insert(
                aggregated_function_matches.end(),
                make_move_iterator(thread_list.begin()),
                make_move_iterator(thread_list.end())
            );
        }
    }

    return aggregated_function_matches;
}

vector<vector<int>> report_long_matches_sl_tail(
    int site_index_0_based,
    int n_hap,
    int num_threads,
    int llm_len,
    const vector<int>& old_pda_parr,
    const vector<int>& old_pda_mlens)
{
    // Tối ưu bằng cách dự đoán kích thước
    const int max_expected_matches = min(n_hap * n_hap / 4, 100000);
    vector<vector<vector<int>>> per_thread_matches_storage(num_threads);
    for (auto& thread_storage : per_thread_matches_storage) {
        thread_storage.reserve(max_expected_matches / num_threads);
    }

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(dynamic, 32)
        for (int h = 0; h < n_hap - 1; ++h) {
            int min_val = numeric_limits<int>::max();
            for (int i = h + 1; i < n_hap; ++i) {
                if (old_pda_mlens[old_pda_parr[i]] < llm_len) {
                    break;
                }
                min_val = min(min_val, old_pda_mlens[old_pda_parr[i]]);
                per_thread_matches_storage[thread_id].push_back({
                    site_index_0_based + 1,
                    min(old_pda_parr[h], old_pda_parr[i]),
                    max(old_pda_parr[h], old_pda_parr[i]),
                    min_val
                });
            }
        }
    }

    // Tính toán tổng số matches để dự đoán kích thước
    size_t total_matches = 0;
    for (const auto& thread_list : per_thread_matches_storage) {
        total_matches += thread_list.size();
    }

    vector<vector<int>> aggregated_function_matches;
    aggregated_function_matches.reserve(total_matches);

    // Sử dụng move iterator để tối ưu việc nối
    for (auto& thread_list : per_thread_matches_storage) {
        if (!thread_list.empty()) {
            aggregated_function_matches.insert(
                aggregated_function_matches.end(),
                make_move_iterator(thread_list.begin()),
                make_move_iterator(thread_list.end())
            );
        }
    }

    return aggregated_function_matches;
}

vector<vector<int>> build_and_match(
    const vector<vector<int>>& haplotypes_sites_matrix,
    int num_threads,
    int L_long_match_min_len)
{
    if (haplotypes_sites_matrix.empty() || haplotypes_sites_matrix[0].empty()) {
        return {};
    }

    const int n_hap = haplotypes_sites_matrix.size();
    const int n_sites = haplotypes_sites_matrix[0].size();

    int block_size = (n_hap + num_threads - 1) / num_threads;
    if (num_threads <= 0) block_size = n_hap;

    // Khởi tạo các vector với kích thước cố định
    vector<int> new_pda_parr(n_hap);
    vector<int> new_pda_mlens(n_hap);
    int new_pda_zerocnt = 0;

    vector<int> old_pda_parr(n_hap);
    vector<int> old_pda_mlens(n_hap);
    int old_pda_zerocnt = 0;

    vector<int> ps_holder(n_hap);
    vector<int> offsets_holder(max(1, num_threads));

    // Dự đoán kích thước cho all_matches_report
    vector<vector<int>> all_matches_report;
    all_matches_report.reserve(min(n_hap * n_sites / 10, 1000000));

    // Khởi tạo sắp xếp ban đầu cho site đầu tiên
    initial_sort(haplotypes_sites_matrix, 0, n_hap, num_threads, block_size,
                old_pda_parr, old_pda_mlens, old_pda_zerocnt,
                ps_holder, offsets_holder);

    for (int site_idx = 1; site_idx < n_sites; ++site_idx) {
        // Xây dựng mảng p mới
        core_p_arr(haplotypes_sites_matrix, site_idx, n_hap, num_threads, block_size,
                  old_pda_parr, new_pda_parr, new_pda_zerocnt,
                  ps_holder, offsets_holder);

        // Xây dựng mảng d mới
        core_d_arr(haplotypes_sites_matrix, site_idx, n_hap, num_threads, block_size,
                  old_pda_parr, old_pda_mlens, new_pda_mlens, ps_holder);

        // Báo cáo các chuỗi khớp dài
        vector<vector<int>> current_site_matches = report_long_matches_sl(
            site_idx - 1, haplotypes_sites_matrix, site_idx, n_hap, num_threads,
            L_long_match_min_len, old_pda_parr, old_pda_mlens);

        // Nối kết quả
        if (!current_site_matches.empty()) {
            all_matches_report.insert(
                all_matches_report.end(),
                make_move_iterator(current_site_matches.begin()),
                make_move_iterator(current_site_matches.end())
            );
        }

        // Hoán đổi các mảng cũ và mới
        old_pda_parr.swap(new_pda_parr);
        old_pda_mlens.swap(new_pda_mlens);
        old_pda_zerocnt = new_pda_zerocnt;
    }

    // Xử lý site cuối cùng
    if (n_sites > 0) {
        vector<vector<int>> tail_matches = report_long_matches_sl_tail(
            n_sites - 1, n_hap, num_threads, L_long_match_min_len,
            old_pda_parr, old_pda_mlens);

        // Nối kết quả
        if (!tail_matches.empty()) {
            all_matches_report.insert(
                all_matches_report.end(),
                make_move_iterator(tail_matches.begin()),
                make_move_iterator(tail_matches.end())
            );
        }
    }

    return all_matches_report;
}



std::vector<std::vector<int>> build_and_match_maximal(
    const std::vector<std::vector<int>>& hap_map_original,
    int num_threads,
    int L_min_len
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

    // Create transposed matrix for better memory access patterns (parallelized)
    std::vector<std::vector<int>> Xt(N, std::vector<int>(M));
    #pragma omp parallel for num_threads(num_threads)
    for (int c = 0; c < N; ++c) {
        for (int r = 0; r < M; ++r) {
            Xt[c][r] = hap_map_original[r][c];
        }
    }

    // Process sites in parallel chunks
    const int chunk_size = std::max(1, N / (num_threads * 2)); // Adjust chunk size based on threads
    const int num_chunks = (N + chunk_size - 1) / chunk_size;
    
    // Each thread will process its own chunk of sites and produce partial results
    std::vector<std::vector<std::vector<int>>> thread_matches(num_chunks);
    
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (int chunk = 0; chunk < num_chunks; ++chunk) {
        int start_k = chunk * chunk_size;
        int end_k = std::min(N, (chunk + 1) * chunk_size);
        
        // Each thread has its own copy of PBWT arrays
        std::vector<int> local_a(M);
        std::iota(local_a.begin(), local_a.end(), 0);
        std::vector<int> local_b(M), local_d(M, 0), local_e(M);
        
        // Initialize PBWT arrays for this chunk
        // If not the first chunk, we need to compute the correct starting state
        if (start_k > 0) {
            // Sequential build up to the starting point of this chunk
            for (int k = 0; k < start_k; ++k) {
                algorithm_2_BuildPrefixAndDivergenceArrays(Xt[k], k, local_a, local_b, local_d, local_e, M);
            }
        }
        
        // Process sites in this chunk
        std::vector<std::vector<int>> chunk_matches;
        for (int k = start_k; k < end_k; ++k) {
            auto k_matches = algorithm_4_ReportMaximalMatches(hap_map_original, local_a, local_d, k);
            chunk_matches.insert(chunk_matches.end(), k_matches.begin(), k_matches.end());
            
            if (k < N - 1) {
                algorithm_2_BuildPrefixAndDivergenceArrays(Xt[k], k, local_a, local_b, local_d, local_e, M);
            }
        }
        
        thread_matches[chunk] = std::move(chunk_matches);
    }
    
    // Combine results from all threads
    std::vector<std::vector<int>> maximal_matches;
    const int estimated_initial_matches = std::min(M * N / 10, 10000);
    maximal_matches.reserve(estimated_initial_matches);
    
    for (auto& chunk_matches : thread_matches) {
        maximal_matches.insert(maximal_matches.end(), chunk_matches.begin(), chunk_matches.end());
    }

    // Sort and remove duplicates - exactly as in pbwt-single.cpp
    if (!maximal_matches.empty()) {
        std::sort(maximal_matches.begin(), maximal_matches.end());
        maximal_matches.erase(std::unique(maximal_matches.begin(), maximal_matches.end()), maximal_matches.end());
    }

    return maximal_matches;
}

// Implementation of Algorithm 4 from the paper
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

int main(int argc, char *argv[]) {
    parallel_run(argc, argv);
}