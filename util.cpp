#include "util.h"
#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <fstream>
#include <iterator>
#include <bits/stdc++.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;
using namespace boost::iostreams;
using namespace chrono;
// Khai báo nguyên mẫu các hàm từ pbwt-single.cpp

void algorithm_2_BuildPrefixAndDivergenceArrays(
    const std::vector<int>& x_k,
    int k,
    std::vector<int>& a,
    std::vector<int>& b,
    std::vector<int>& d,
    std::vector<int>& e,
    int M_haps
);

std::vector<std::vector<int>> algorithm_4_ReportMaximalMatches(
    const std::vector<std::vector<int>>& hap_matrix,
    const std::vector<int>& p_k,
    const std::vector<int>& d_k,
    int k
);

void parallel_run(int argc, char *argv[]) {
    string f = argv[1];
   // printf("Read file: %s\n", f.c_str());
    vector<vector<int> > X = read_hap(f);

    const int retry = 1;
    const int L = 4;
   // printf("Test size: %lu haps x %lu sites\n", X.size(), X[0].size());
    //printf("Retry: %d\n", retry);
   // printf("Matched length: %d\n", L);

    //int nthreads = static_cast<int>(thread::hardware_concurrency());
    int nthreads = 4;
    //printf("Max thread:%d\n\n", nthreads);
    map<tuple<int, int, int>, vector<int>> match_map;
    vector<vector<int> > matches_l;
    printf("========== Long Matches ==========\n");
    for (int i = 2; i <= nthreads; i++) {
        // printf("Num thread:%d\n", i);
       long long duration_total_run[retry];
        //signed long int duration_total_run = 0;
        high_resolution_clock::time_point start, stop;
        vector<vector<int> >matches_l = build_and_match(X, i, L);
        vector<vector<vector<int> > > res;
        for (int j = 0; j < retry; ++j) {
            start = high_resolution_clock::now();
             matches_l = build_and_match(X, i, L);
            //  matches = build_and_match_maximal(X, i, 1);
            stop = high_resolution_clock::now();
            duration_total_run[j-1] = duration_cast<microseconds>(stop - start).count();
            //printf("Time: %ld us\n", duration_total_run[j-1]);
        }
        // for (const auto& m : matches) {
        //     auto key = make_tuple(m[0], m[1], m[2]); // (pos, start, end)
        //     auto it = match_map.find(key);
        //     if (it == match_map.end() || m[3] > it->second[3]) {
        //         match_map[key] = m; // Giữ đoạn có length lớn nhất
        //     }
        // }
        
        
        print_matches(matches_l);
        // sort duration_total_run
        sort(duration_total_run, duration_total_run+retry);

    }
    // vector<vector<int>> merged_matches;
    // for (const auto& kv : match_map) {
    //     merged_matches.push_back(kv.second);
    // }

    // sort(merged_matches.begin(), merged_matches.end(), [](const vector<int>& a, const vector<int>& b) {
    //     return a[0] == b[0] ? (a[1] == b[1] ? a[2] < b[2] : a[1] < b[1]) : a[0] < b[0];
    // });

    
    vector<vector<int> >matches_max = build_and_match_maximal(X, 5, 1);
    printf("\n========== Maximal Matches ==========\n");
    print_matches(matches_max);
}

void single_run(int argc, char *argv[]) {
    string f = argv[1];
    vector<vector<int>> X = read_hap(f);

    const int retry = 5;
    const int L = 1;
    int nthreads = static_cast<int>(thread::hardware_concurrency());

    vector<long long> duration_total_run;
    high_resolution_clock::time_point start, stop;

    // ===================== PBWT LLM =======================
    vector<vector<int>> llm_matches, maximal_matches;

    for (int j = 0; j < retry; ++j) {
        start = high_resolution_clock::now();
        llm_matches = build_and_match(X, 1, 4);  // Long matches
        stop = high_resolution_clock::now();
        duration_total_run.push_back(duration_cast<microseconds>(stop - start).count());
    }

    printf("========== Long Matches ==========\n");
    print_matches(llm_matches);

    printf("\n========== Maximal Matches ==========\n");
    auto maximal = build_and_match_maximal(X, 1, 1);
    print_matches(maximal);

}


void print_matches(const vector<vector<int> > &all_matches_report)
{
   vector<vector<int> > sorted =(all_matches_report);

    // Sort the filtered matches by column 1 (position)
    sort(sorted.begin(), sorted.end(),
         [](const vector<int>& a, const vector<int>& b) {
             return a[0] == b[0] ? a[1] == b[1] ? a[2] < b[2] : a[1] < b[1] : a[0] < b[0];
         });
   
    // Print the sorted matches
    for (const auto &match: sorted) {
        printf("Position %d: [%d, %d]\n",
               match[0], match[1], match[2]);

             
    }
    
}

vector<vector<int> > read_hap(const string &filename) {
    if (filename.find(".gz") != string::npos)
        return read_hap_gz(filename);
    return read_hap_normal(filename);
}

vector<vector<int> > read_hap_normal(const string &filename) {
    vector<vector<int> > res;
    ifstream ifs(filename);       // open the file
    string tempstr;
    while (getline(ifs, tempstr)) {
        stringstream lineStream(tempstr);
        vector<int> numbers(istream_iterator<int>(lineStream),{});
        if (res.empty()) {
            for (auto number: numbers) {
                res.push_back(vector<int>{number});
            }
        }
        else {
            for (int i = 0; i < numbers.size(); ++i) {
                res[i].push_back(numbers[i]);
            }
        }
    }
    ifs.close();
    return res;
}

vector<vector<int> > read_hap_gz(const string &filename){
    vector<vector<int> > res;
    ifstream file(filename, ios_base::in | ios_base::binary);
    filtering_streambuf<input> inbuf;
    inbuf.push(gzip_decompressor());
    inbuf.push(file);
    //Convert streambuf to istream
    istream instream(&inbuf);
    //Iterate lines
    string line;
    while(getline(instream, line)) {
        stringstream lineStream(line);
        vector<int> numbers(istream_iterator<int>(lineStream), {});
        if (res.empty()) {
            for (auto number: numbers) {
                res.push_back(vector<int>{number});
            }
        } else {
            for (int i = 0; i < numbers.size(); ++i) {
                res[i].push_back(numbers[i]);
            }
        }
    }
    file.close();
    return res;
}
