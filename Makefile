build:
	g++ pbwt.cpp -o pbwt -fopenmp -lz -lboost_iostreams -std=c++11 -I"${CONDA_PREFIX}"/include
	g++ randomizer.cpp -o randomizer -fopenmp -std=c++11 -I"${CONDA_PREFIX}"/include
	g++ pbwt_1.cpp -o pbwt_1 -lz -lboost_iostreams -std=c++11 -I"${CONDA_PREFIX}"/include
