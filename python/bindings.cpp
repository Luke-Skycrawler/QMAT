#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <string>
#include "pyqmat.h"
#include "point_adder.h"
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

    py::class_<SlabMesh>(m, "SlabMesh");
    
    py::class_<QMAT, SlabMesh>(m, "qmat")
        .def(py::init<const std::string &, const std::string &>())
        .def(py::init<const QMAT &>())
        .def(py::init<const QMAT &, const std::string &>())
        .def("hausdorff", &QMAT::hausdorff)
        .def("simplify_slab", &QMAT::simplify_slab)
        .def("export_ply", &QMAT::export_ply)
        .def("init_mergelist", &QMAT::initMergeList)
        .def("init_collapse_queue", &QMAT::initCollapseQueue)
        .def("simplify", &QMAT::Simplify)
        .def("adjust_storage", &QMAT::AdjustStorage)
        .def("export_mergelist", &QMAT::exportMergeList)
        .def("export_ma", &QMAT::export_ma)
        .def("get_diagonal", &QMAT::get_diagonal);

    py::class_<PointAdder>(m, "PointAdder")
        .def(py::init<double, SlabMesh &, SlabMesh &>())
        .def("set_mergelist", &PointAdder::set_collapsed_list)
        .def("add_new_node", &PointAdder::add_new_noded)
        .def("export_ply", &PointAdder::export_ply);
}