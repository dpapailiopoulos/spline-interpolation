#ifndef READ_PATH_H
#define READ_PATH_H

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <eigen3/Eigen/Core>


Eigen::MatrixXd readTrack(const std::string& filename) {
    int cols{ 0 };
    int rows{ 0 };
    int i{ 0 };
    std::array<double, (int)1e6> buffer{};
    std::string line;

    // Read numbers from the file into the buffer
    std::ifstream data_file;
    data_file.open(filename);
    while (!data_file.eof())
    {
        // Get next line and if empty stop iterating through the file
        std::getline(data_file, line);
        if (line == "") break;
        
        // Convert line to stream to extract point values
        std::stringstream stream(line);
        while (!stream.eof()) {
            stream >> buffer[i];
            i++;
        }

        line.clear();
        // Set the number of columns at the first iteration
        if (cols == 0) cols = i;
        rows++;
    }
    data_file.close();

    // Populate the eigen matrix with the numbers from the buffer
    Eigen::MatrixXd result{ rows, cols };
    for (int i{ 0 }; i < rows; i++)
        for (int j{ 0 }; j < cols; j++)
            result(i, j) = buffer[cols * i + j];
    return result;
}
#endif