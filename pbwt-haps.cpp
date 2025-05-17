#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include <omp.h>
#include <cmath>
#include <limits>
#include <string>

#include "util.h"

using namespace std;
using namespace chrono;

void core_p_arr(
    int site_idx, // Not used in C# version's coreP_Arr for map access, map passed directly
    const function<char(int)>& get_hap_val, // Corresponds to map.Get_HapVal(hap_id)
    int n_hap, // Program.nHap
    int num_threads, // Program.nThread
    int block_size,
    const vector<int>& old_pda_parr, // oldPDA.pArr
    vector<int>& new_pda_parr, // newPDA.pArr
    int& new_pda_zerocnt, // newPDA.zeroCnt
    vector<int>& ps_holder, // psHolder
    vector<int>& offsets_holder) // offsetsHolder
{
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 0; i < num_threads; ++i) {
        int s = i * block_size;
        if (s >= n_hap) continue;
        int e = min(s + block_size, n_hap);

        if (s < e) {
            if (get_hap_val(old_pda_parr[s]) == '0') {
                ps_holder[s] = 0;
            } else {
                ps_holder[s] = 1;
            }

            for (int k = s + 1; k < e; ++k) {
                if (get_hap_val(old_pda_parr[k]) == '0') {
                    ps_holder[k] = ps_holder[k - 1];
                } else {
                    ps_holder[k] = ps_holder[k - 1] + 1;
                }
            }
            offsets_holder[i] = ps_holder[e - 1];
        } else {
             offsets_holder[i] = 0;
        }
    }

    for (int i = 1; i < num_threads; ++i) {
        offsets_holder[i] += offsets_holder[i - 1];
    }

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 1; i < num_threads; ++i) {
        int s = i * block_size;
        if (s >= n_hap) continue;
        int e = min(s + block_size, n_hap);
        for (int k = s; k < e; ++k) {
            ps_holder[k] += offsets_holder[i - 1];
        }
    }

    if (n_hap > 0) {
        new_pda_zerocnt = n_hap - ps_holder[n_hap - 1];
    } else {
        new_pda_zerocnt = 0;
    }
    int one_off = new_pda_zerocnt -1; // C# oneOff

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 0; i < n_hap; ++i) {
        if (get_hap_val(old_pda_parr[i]) == '0') {
            new_pda_parr[i - ps_holder[i]] = old_pda_parr[i];
        } else {
            new_pda_parr[ps_holder[i] + one_off] = old_pda_parr[i];
        }
    }
}

void initial_sort(
    const function<char(int)>& get_hap_val_initial, // map.Get_HapVal(hap_id) for initial site
    int n_hap, // Program.nHap
    int num_threads, // Program.nThread
    int block_size,
    vector<int>& old_pda_parr, // oldPDA.pArr (output)
    vector<int>& old_pda_mlens, // oldPDA.mLens (output)
    int& old_pda_zerocnt, // oldPDA.zeroCnt (output)
    vector<int>& ps_holder, // psHolder
    vector<int>& offsets_holder) // offsetsHolder
{
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 0; i < num_threads; ++i) {
        int s = i * block_size;
        if (s >= n_hap) continue;
        int e = min(s + block_size, n_hap);

        if (s < e) {
            if (get_hap_val_initial(s) == '0') { // Accessing by original hap_id 's'
                ps_holder[s] = 0;
            } else {
                ps_holder[s] = 1;
            }

            for (int k = s + 1; k < e; ++k) { // Accessing by original hap_id 'k'
                if (get_hap_val_initial(k) == '0') {
                    ps_holder[k] = ps_holder[k - 1];
                } else {
                    ps_holder[k] = ps_holder[k - 1] + 1;
                }
            }
            offsets_holder[i] = ps_holder[e - 1];
        } else {
            offsets_holder[i] = 0;
        }
    }

    for (int i = 1; i < num_threads; ++i) {
        offsets_holder[i] += offsets_holder[i - 1];
    }

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 1; i < num_threads; ++i) {
        int s = i * block_size;
        if (s >= n_hap) continue;
        int e = min(s + block_size, n_hap);
        for (int k = s; k < e; ++k) {
            ps_holder[k] += offsets_holder[i - 1];
        }
    }

    if (n_hap > 0) {
        old_pda_zerocnt = n_hap - ps_holder[n_hap - 1];
    } else {
        old_pda_zerocnt = 0;
    }
    int one_off = old_pda_zerocnt -1;

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 0; i < n_hap; ++i) { // i is original haplotype ID
        if (get_hap_val_initial(i) == '0') {
            old_pda_parr[i - ps_holder[i]] = i;
            // old_pda_mlens[i - ps_holder[i]] = 1; // C# original: oldPDA.mLens[i - psHolder[i]] = 1;
                                                 // This should be oldPDA.mLens[original_hap_id]
            old_pda_mlens[i] = 1; // Corrected: D-value is for the haplotype ID itself
        } else {
            old_pda_parr[ps_holder[i] + one_off] = i;
            // old_pda_mlens[ps_holder[i] + one_off] = 1;
            old_pda_mlens[i] = 1; // Corrected
        }
    }

    // D Arr adjust
    if (n_hap > 0) { // Ensure pArr is not empty
        if (old_pda_parr.size() > 0) old_pda_mlens[old_pda_parr[0]] = 0; // C# oldPDA.mLens[oldPDA.pArr.First()] = 0;
        if (old_pda_zerocnt != n_hap) {
            if (static_cast<int>(old_pda_parr.size()) > one_off + 1 && one_off + 1 >= 0) { // Check bounds
                old_pda_mlens[old_pda_parr[one_off + 1]] = 0; // C# oldPDA.mLens[oldPDA.pArr[oneOff + 1]] = 0;
            }
        }
    }
}

void core_d_arr(
    const function<char(int)>& get_hap_val, // map.Get_HapVal
    int n_hap, // Program.nHap
    int num_threads, // Program.nThread
    int block_size,
    const vector<int>& old_pda_parr, // oldPDA.pArr
    const vector<int>& old_pda_mlens, // oldPDA.mLens
    vector<int>& new_pda_mlens, // newPDA.mLens
    const vector<int>& ps_holder // psHolder (used in C# for block analysis, simplified here)
)
{
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int i = 0; i < num_threads; ++i) { // i is block_idx
        int s = i * block_size;
        if (s >= n_hap) continue;
        int e = min(s + block_size, n_hap);

        int prv_low_m_zero = -1;
        int prv_low_m_one = -1;
        int h_id;
        // int seek_hid; // Not directly used in this simplified C++ translation

        if (s == 0) { // First block
            // prv_low_m_zero = -1; // Already initialized
            // prv_low_m_one = -1;  // Already initialized
            for (int k = s; k < e; ++k) {
                h_id = old_pda_parr[k];
                if (get_hap_val(h_id) == '0') {
                    new_pda_mlens[h_id] = min(old_pda_mlens[h_id], prv_low_m_one) + 1;
                    prv_low_m_one = numeric_limits<int>::max(); // C# int.MaxValue
                    prv_low_m_zero = min(prv_low_m_zero, old_pda_mlens[h_id]);
                } else { // incoming is 1
                    new_pda_mlens[h_id] = min(old_pda_mlens[h_id], prv_low_m_zero) + 1;
                    prv_low_m_zero = numeric_limits<int>::max(); // C# int.MaxValue
                    prv_low_m_one = min(prv_low_m_one, old_pda_mlens[h_id]);
                }
            }
        } else { // Other blocks
            // Simplified lookup from C#; a full translation of the prvLowM search is complex
            // and requires careful handling of psHolder or similar logic.
            // For this version, we'll do a simpler (less efficient) backward scan.
            prv_low_m_zero = -1;
            prv_low_m_one = -1;

            bool min_zero_search = true;
            bool min_one_search = true;

            if (ps_holder[s - 1] == s) { // A there is no 0 in upper blocks
                prv_low_m_zero = -1;
                min_zero_search = false;
            }
            if (ps_holder[s - 1] == 0) { // B there is no 1 in upper blocks
                prv_low_m_one = -1;
                min_one_search = false;
            }

            // This is a simplification of the C# prvLowM search
            if (min_zero_search) {
                int seek_index = s - 1;
                while(seek_index >= 0 && get_hap_val(old_pda_parr[seek_index]) != '0') {
                    seek_index--;
                }
                if (seek_index >= 0) {
                    prv_low_m_zero = old_pda_mlens[old_pda_parr[seek_index]];
                     while (seek_index >= 0 && get_hap_val(old_pda_parr[seek_index]) == '0') {
                        prv_low_m_zero = min(old_pda_mlens[old_pda_parr[seek_index]], prv_low_m_zero);
                        seek_index--;
                    }
                    // C# has an additional check: if (seekIndex == -1) prvLowM_Zero = -1;
                    // This case is implicitly handled if the loop doesn't find any '0'.
                } else {
                    prv_low_m_zero = -1; // No '0' found above
                }
            }

            if (min_one_search) {
                int seek_index = s - 1;
                while(seek_index >= 0 && get_hap_val(old_pda_parr[seek_index]) == '0') { // C# map.Get_HapVal(oldPDA.pArr[seekIndex]) == '0'
                    seek_index--;
                }
                if (seek_index >=0) {
                    prv_low_m_one = old_pda_mlens[old_pda_parr[seek_index]];
                    while (seek_index >= 0 && get_hap_val(old_pda_parr[seek_index]) != '0') { // C# map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0'
                        prv_low_m_one = min(old_pda_mlens[old_pda_parr[seek_index]], prv_low_m_one);
                        seek_index--;
                    }
                } else {
                    prv_low_m_one = -1; // No '1' found above
                }
            }


            for (int k = s; k < e; ++k) {
                h_id = old_pda_parr[k];
                if (get_hap_val(h_id) == '0') {
                    new_pda_mlens[h_id] = min(old_pda_mlens[h_id], prv_low_m_one) + 1;
                    prv_low_m_one = numeric_limits<int>::max();
                    prv_low_m_zero = min(prv_low_m_zero, old_pda_mlens[h_id]);
                } else { // incoming is 1
                    new_pda_mlens[h_id] = min(old_pda_mlens[h_id], prv_low_m_zero) + 1;
                    prv_low_m_zero = numeric_limits<int>::max();
                    prv_low_m_one = min(prv_low_m_one, old_pda_mlens[h_id]);
                }
            }
        }
    }
}


void report_long_matches_sl( // Renamed to match C#
    int site_index, // siteIndex
    const function<char(int)>& get_hap_val, // map.Get_HapVal
    int n_hap, // Program.nHap
    int num_threads, // Program.nThread (used for omp_set_num_threads)
    int llm_len, // Program.LLM_Len
    const vector<int>& old_pda_parr, // oldPDA.pArr
    const vector<int>& old_pda_mlens, // oldPDA.mLens
    vector<vector<int>>& all_matches_report) // Output like Program.BW.Add
{
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for schedule(dynamic)
    for (int h = 0; h < n_hap - 1; ++h) {
        int min_val = numeric_limits<int>::max();
        for (int i = h + 1; i < n_hap; ++i) {
            if (old_pda_mlens[old_pda_parr[i]] < llm_len) { // C# oldPDA.mLens[oldPDA.pArr[i]]
                break;
            }
            min_val = min(min_val, old_pda_mlens[old_pda_parr[i]]);

            if (get_hap_val(old_pda_parr[h]) != get_hap_val(old_pda_parr[i])) {
                 #pragma omp critical
                 {
                    // C# Program.BW.Add(oldPDA.pArr[h] + "\t" + oldPDA.pArr[i] + "\t" + siteIndex + "\t" + minVal);
                    // Output format: {hap1, hap2, siteIndex, minVal}
                    // The problem asks for {site, hap1, hap2} - minVal is implicitly the L in L-long.
                    // If minVal is the actual length, it should be included.
                    // Let's assume the request means {site, hap1, hap2, length}
                    all_matches_report.push_back({min(old_pda_parr[h], old_pda_parr[i]),
                                                   max(old_pda_parr[h], old_pda_parr[i]),
                                                   site_index,
                                                   min_val});
                 }
            }
        }
    }
}

void report_long_matches_sl_tail( // Renamed to match C#
    int site_index, // siteIndex
    int n_hap, // Program.nHap
    int num_threads,
    int llm_len, // Program.LLM_Len
    const vector<int>& old_pda_parr, // oldPDA.pArr
    const vector<int>& old_pda_mlens, // oldPDA.mLens
    vector<vector<int>>& all_matches_report)
{
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for schedule(dynamic)
    for (int h = 0; h < n_hap - 1; ++h) {
        int min_val = numeric_limits<int>::max();
        for (int i = h + 1; i < n_hap; ++i) {
            if (old_pda_mlens[old_pda_parr[i]] < llm_len) {
                 break;
            }
            min_val = min(min_val, old_pda_mlens[old_pda_parr[i]]);
            #pragma omp critical
            {
                all_matches_report.push_back({min(old_pda_parr[h], old_pda_parr[i]),
                                               max(old_pda_parr[h], old_pda_parr[i]),
                                               site_index,
                                               min_val});
            }
        }
    }
}


vector<vector<int>> build_and_match(
    const vector<vector<int>>& haplotypes_sites_matrix,
    int num_threads,
    int L_long_match_min_len) // Program.LLM_Len
{
    if (haplotypes_sites_matrix.empty() || haplotypes_sites_matrix[0].empty()) {
        return {};
    }
    int n_hap = haplotypes_sites_matrix.size(); // Program.nHap
    int n_sites = haplotypes_sites_matrix[0].size();

    // omp_set_num_threads(num_threads); // Set inside each parallel region
    int block_size = (n_hap + num_threads - 1) / num_threads;
    if (num_threads <= 0) block_size = n_hap;

    vector<int> new_pda_parr(n_hap); // newPDA.pArr
    vector<int> new_pda_mlens(n_hap); // newPDA.mLens
    int new_pda_zerocnt = 0; // newPDA.zeroCnt

    vector<int> old_pda_parr(n_hap); // oldPDA.pArr
    vector<int> old_pda_mlens(n_hap); // oldPDA.mLens
    int old_pda_zerocnt = 0; // oldPDA.zeroCnt

    vector<int> ps_holder(n_hap); // psHolder
    vector<int> offsets_holder(max(1,num_threads)); // offsetsHolder

    vector<vector<int>> all_matches_report; // Output like Program.BW

    // Lambda to get allele for current site being processed by PBWT
    // site_idx is the current SNP/locus index in haplotypes_sites_matrix
    auto get_hap_val_for_site_k =
        [&](int current_processing_site_idx, int original_hap_id) -> char {
        return haplotypes_sites_matrix[original_hap_id][current_processing_site_idx] == 0 ? '0' : '1';
    };

    // InitialSort uses the alleles at the first site (site 0)
    initial_sort(
        [&](int original_hap_id){ return get_hap_val_for_site_k(0, original_hap_id); },
        n_hap, num_threads, block_size,
        old_pda_parr, old_pda_mlens, old_pda_zerocnt,
        ps_holder, offsets_holder);

    for (int site_idx = 1; site_idx < n_sites; ++site_idx) { // site_idx is the current site for PBWT update
        // coreP_Arr
        core_p_arr(
            site_idx, // Not strictly needed if lambda captures it
            [&](int original_hap_id){ return get_hap_val_for_site_k(site_idx, original_hap_id); },
            n_hap, num_threads, block_size,
            old_pda_parr,
            new_pda_parr, new_pda_zerocnt,
            ps_holder, offsets_holder);

        // coreD_Arr
        // ps_holder from the core_p_arr immediately above is relevant for complex D array calculation
        // The C# code uses psHolder which reflects the '0'/'1' counts *after* sorting by site_idx.
        // So, the ps_holder from the most recent core_p_arr is needed for a more accurate core_d_arr.
        core_d_arr(
            [&](int original_hap_id){ return get_hap_val_for_site_k(site_idx, original_hap_id); },
            n_hap, num_threads, block_size,
            old_pda_parr, // pArr from previous site (k-1)
            old_pda_mlens, // mLens from previous site (k-1)
            new_pda_mlens, // mLens for current site (k)
            ps_holder // psHolder from current site's P array computation
        );

        // ReportLongMatches_SL
        // Reports matches ending at site_idx-1, using alleles at site_idx for mismatch check
        // It uses oldPDA from C#, which means p_arr and m_lens *before* the current site's update.
        // In our loop: old_pda_parr and old_pda_mlens are from site_idx-1.
        // new_pda_parr and new_pda_mlens are for site_idx.
        // The C# ReportLongMatches_SL uses `oldPDA`, which means the state *before* OneSort.
        // OneSort does: coreP, coreD, then swaps oldPDA and newPDA.
        // So, for reporting, we should use the state corresponding to *before* the current site_idx update.
        // However, the mismatch check uses `map.Get_HapVal` which is for the *current* site.

        report_long_matches_sl(
            site_idx - 1, // reported_site_idx is the site where the match ends
            [&](int original_hap_id){ return get_hap_val_for_site_k(site_idx, original_hap_id); }, // Alleles at current site_idx for mismatch
            n_hap, num_threads, L_long_match_min_len,
            old_pda_parr, // P array from site_idx-1
            old_pda_mlens, // D array from site_idx-1
            all_matches_report);

        // Swap for next iteration (temPDA = oldPDA; oldPDA = newPDA; newPDA = temPDA;)
        old_pda_parr.swap(new_pda_parr);
        old_pda_mlens.swap(new_pda_mlens);
        old_pda_zerocnt = new_pda_zerocnt;
    }

    if (n_sites > 0) {
        // ReportLongMatches_SL_Tail
        // Uses the state of oldPDA after the last OneSort (i.e., for the last site)
        report_long_matches_sl_tail(
            n_sites - 1,
            n_hap, num_threads, L_long_match_min_len,
            old_pda_parr, // P and D arrays for the last site
            old_pda_mlens,
            all_matches_report);
    }

    // The problem output {site, hap1, hap2}. Current output includes length.
    // Adjusting output format:
    vector<vector<int>> final_report;
    for(const auto& match_with_len : all_matches_report) {
        if (match_with_len.size() == 4) { // site, hap1, hap2, length
            final_report.push_back({match_with_len[2]+1, match_with_len[0], match_with_len[1]});
        }
    }

    return final_report;
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

    int nthreads = static_cast<int>(thread::hardware_concurrency());
    printf("Max thread:%d\n\n", nthreads);
    for (int i = 1; i <= nthreads; ++i) {
        printf("Num thread:%d\n", i);
        signed long int duration_total_build = 0;
        signed long int duration_total_match = 0;
        high_resolution_clock::time_point start, stop;
        vector<vector<int> > matches;
        vector<vector<vector<int> > > res;
        for (int j = 0; j < retry; ++j) {
            start = high_resolution_clock::now();
            matches = build_and_match(X,1,4);
            stop = high_resolution_clock::now();
            duration_total_build += duration_cast<microseconds>(stop - start).count();
            print_matches(matches);
        }
        printf("Total time: %ld us\n\n", duration_total_match / retry + duration_total_build / retry);
    }
}