#ifndef FILEHANDLER_H
#define FILEHANDLER_H

#include <fstream>
#include <stdlib.h>
#include <sstream>
void mkdir(std::string s){
	std::stringstream aux;
	aux << "mkdir -p " + s + "/";
	system(aux.str().c_str());	
}

void storeXYData(double X, double Y, std::string fileName, bool append = true, std::string delimiter = "\t") {
    std::ofstream file;
    if (append) {
        file.open(fileName, std::ios::app);
    } else {
        file.open(fileName);
    }
    file << X << delimiter << Y << std::endl;
    file.close();
}

template <typename T>
void storeData(T data, std::string fileName) {
    std::ofstream file;
    file.open(fileName);
    file << data << std::endl;
    file.close();
}
#endif
