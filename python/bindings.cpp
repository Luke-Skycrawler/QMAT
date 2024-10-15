#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <string>
#include "pyqmat.h"
namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(pyqmat, m)
{
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: shaysweep

        .. autosummary::
           :toctree: _generate


    )pbdoc";

    py::class_<QMAT>(m, "qmat")
        .def(py::init<const std::string&, const std::string&>())
        .def("compute_hausdorff_distance", &QMAT::ComputeHausdorffDistance)
        .def("simplify_slab", &QMAT::simplifySlab)
        .def("export_ply", &QMAT::ExportPly)
        .def("export_ma", &QMAT::ExportMA)
        .def("export_hausdorff_distance", &QMAT::export_hausdorff_distance);
}