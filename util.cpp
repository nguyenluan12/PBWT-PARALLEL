#include "util.h"
#include <vector>
#include <fstream>
#include <iterator>
#include <bits/stdc++.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;
using namespace boost::iostreams;

void print_matches(const vector<vector<int>>& matches) {
    // Copy & sort theo phần tử đầu tiên
    vector<vector<int>> sorted_matches(matches);
    sort(sorted_matches.begin(), sorted_matches.end(),
         [](const vector<int>& a, const vector<int>& b) {
             return a[0] < b[0];
         });

    for (const auto &match: sorted_matches) {
        // cout << "Position " << match[0] << ": ma = [";
        // size_t idx = 1;
        // while (idx < match.size() && match[idx] != -1) {
        //     cout << match[idx++];
        //     if (idx < match.size() && match[idx] != -1) {
        //         cout << ", ";
        //     }
        // }
        // cout << "], mb = [";
        // while (idx < match.size() && match[idx] == -1) {
        //     ++idx; // Skip the -1 separator
        // }
        // while (idx < match.size()) {
        //     cout << match[idx++];
        //     if (idx < match.size()) {
        //         cout << ", ";
        //     }
        // }
        // cout << "]" << endl;
        vector<int> values;
        for (size_t i = 1; i < match.size(); ++i) {
            if (match[i] != -1) {
                values.push_back(match[i]);
            }
        }
        sort(values.begin(), values.end());
        printf("Position %d: [", match[0]);
        for (size_t j = 0; j < values.size(); ++j) {
            if (j) printf(", ");
            printf("%d", values[j]+1);
        }
        printf("]\n");
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