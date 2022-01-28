#include <pybind11/pybind11.h>
//#include <pybind11/embed.h>

#include "specToChrom.h"

namespace py = pybind11;

PYBIND11_MODULE(mzMLTrace, m) {
  py::class_<mzMLTrace::specToChrom>(m, "specToChrom")
    .def(py::init<>())
    .def("set_filename", &mzMLTrace::specToChrom::set_filename)
    .def("set_minimum_intensity", &mzMLTrace::specToChrom::set_minimum_intensity)
    .def("print_filename", &mzMLTrace::specToChrom::print_filename)
    .def("is_set", &mzMLTrace::specToChrom::is_set)
    .def("reset", &mzMLTrace::specToChrom::reset)
    .def("readSpectra", &mzMLTrace::specToChrom::readSpectra)
    .def("findChromatograms", &mzMLTrace::specToChrom::findChromatograms)
    .def("writeChromatograms", &mzMLTrace::specToChrom::writeChromatograms)
    .def("get_chrom", &mzMLTrace::specToChrom::get_chrom)
    .def("get_nchrom", &mzMLTrace::specToChrom::get_nchrom)
    .def("get_chrom_len", &mzMLTrace::specToChrom::get_chrom_len);
}
