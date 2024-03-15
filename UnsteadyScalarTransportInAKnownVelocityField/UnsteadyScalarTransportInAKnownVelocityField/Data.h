#pragma once
#include <iostream>
#include <vector>
#include <fstream>

class Data {
	void save(std::vector<double>& data, std::string filename) {
		std::ofstream file("./data/" + filename + ".txt");
		for (int i = 0; i < data.size(); i++) {
			file << data[i];
			if (i < data.size() - 1) {
				file << " ";
			}
		}
		file.close();
	}
	void save(std::vector<std::vector<double>>& data, std::string filename) {
		std::ofstream file("./data/" + filename + ".txt");
		for (int i = 0; i < data.size(); i++) {
			for (int j = 0; j < data[i].size(); j++) {
				file << data[i][j];
				if (j < data[i].size() - 1) {
					file << " ";
				}
			}
			if (i < data.size() - 1) {
				file << "\n";
			}
		}
		file.close();
	}
public:
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> t;
	std::vector<std::vector<std::vector<double>>> phi;
	void save(std::string index) {
		save(x, "x");
		save(y, "y");
		save(phi[phi.size() - 1], "phi" + index);
	}
};
