// basic file operations
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s <haps> <sites>\n", argv[0]);
        return 1;
    }

    int hap = stoi(argv[1]);
    int sites = stoi(argv[2]);

    // Generate the random binary matrix and store it in memory
    vector<vector<int>> matrix(sites, vector<int>(hap));

    for (int j = 0; j < sites; ++j) {
        for (int k = 0; k < hap; ++k) {
            matrix[j][k] = rand() % 2;
        }
    }

    // Output the sitesXhaps matrix (sites as rows, haps as columns)
    ofstream sitesXhapsFile;
    sitesXhapsFile.open(to_string(sites) + "SitesX" + to_string(hap) + "Haps.hap", std::ios::out | std::ios::trunc);
    for (int j = 0; j < sites; ++j) {
        for (int k = 0; k < hap; ++k) {
            sitesXhapsFile << matrix[j][k] << " ";
        }
        sitesXhapsFile << "\n";
    }
    sitesXhapsFile.close();

    // Output the hapsXsites matrix (haps as rows, sites as columns)
    ofstream hapsXsitesFile;
    hapsXsitesFile.open(to_string(hap) + "HapsX" + to_string(sites) + "Sites.hap", std::ios::out | std::ios::trunc);
    for (int k = 0; k < hap; ++k) {
        for (int j = 0; j < sites; ++j) {
            hapsXsitesFile << matrix[j][k] << " ";
        }
        hapsXsitesFile << "\n";
    }
    hapsXsitesFile.close();

    printf("Generated %dx%d matrix and saved to two files.\n", sites, hap);
    return 0;
}

