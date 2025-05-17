// basic file operations
#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char *argv[]) {
    ofstream myfile;
    int hap=stoi(argv[1]);
    int sites=stoi(argv[2]);
    myfile.open(to_string(hap)+"HapX"+to_string(sites)+"Site.hap",std::ios::out | std::ios::trunc);
    for (int j = 0; j < sites; ++j) {
        for (int k = 0; k < hap; ++k) {
            myfile << rand() % 2;
            myfile << " ";
        }
        myfile << "\n";
    }
    myfile.close();
}
