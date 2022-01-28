cmake_minimum_required(VERSION 2.8.12)
project(mzMLTrace)

find_package(pybind11 REQUIRED)
pybind11_add_module(mzMLTrace mzMLTrace.cpp specToChrom.cc spectrum.cc base64.cc chromatogram.cc)
target_link_libraries(mzMLTrace PRIVATE z)
