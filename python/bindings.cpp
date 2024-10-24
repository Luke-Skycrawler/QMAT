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
        .def("hausdorff", &QMAT::hausdorff)
        .def("simplify_slab", &QMAT::simplify_slab)
        .def("export_ply", &QMAT::export_ply)
        .def("init_mergelist", &QMAT::initMergeList)
        .def("init_collapse_queue", &QMAT::initCollapseQueue)
        .def("simplify", &QMAT::Simplify)
        .def("adjust_storage", &QMAT::AdjustStorage)
        .def("export_mergelist", &QMAT::exportMergeList)
        .def("export_ma", &QMAT::export_ma);
}