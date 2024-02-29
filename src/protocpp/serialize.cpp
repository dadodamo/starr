//
// Created by Daniel Adamovic on 28/10/23.
//
#include "serialize.h"

void proto::serialize(sampler_data::samples& sample_stream){
    std::string serialized_samples;
    if (!sample_stream.SerializeToString(&serialized_samples)) {
        std::cerr << "Failed to write sampled data." << std::endl;
    }
    std::ofstream outputFile_samples("samples_serialized.bin", std::ios::binary);
    if (outputFile_samples.is_open()) {
        outputFile_samples.write(serialized_samples.c_str(), serialized_samples.size());
        outputFile_samples.close();
        std::cout << "Serialized data written to samples_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }
}

void proto::serialize_y(y_data::full_y& y_stream) {
    std::string serialized_y;
    if (!y_stream.SerializeToString(&serialized_y)) {
        std::cerr << "Failed to write y data." << std::endl;
    }
    std::ofstream outputFile_samples("y_serialized.bin", std::ios::binary);
    if (outputFile_samples.is_open()) {
        outputFile_samples.write(serialized_y.c_str(), serialized_y.size());
        outputFile_samples.close();
        std::cout << "Serialized data written to y_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }
}
