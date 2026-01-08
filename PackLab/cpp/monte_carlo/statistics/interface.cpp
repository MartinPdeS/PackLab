#include "statistics.h"

#include <pybind11/pybind11.h>


PYBIND11_MODULE(interface_statistics, module) {
    pybind11::class_<Statistics>(module, "Statistics")
        .def_readonly("attempted_insertions", &Statistics::attempted_insertions)
        .def_readonly("accepted_insertions", &Statistics::accepted_insertions)
        .def_readonly("rejected_insertions", &Statistics::rejected_insertions)
        .def_readonly("consecutive_rejections", &Statistics::consecutive_rejections)
        .def_readonly("sphere_count", &Statistics::sphere_count)
        .def_readonly("packing_fraction_geometry", &Statistics::packing_fraction_geometry)
        .def_readonly("packing_fraction_simulator", &Statistics::packing_fraction_simulator)
        .def_readonly("radius_min", &Statistics::radius_min)
        .def_readonly("radius_max", &Statistics::radius_max)
        .def_readonly("radius_mean", &Statistics::radius_mean)
        .def_readonly("radius_median", &Statistics::radius_median)
        .def_readonly("radius_std", &Statistics::radius_std)
        .def_readonly("total_runtime_seconds", &Statistics::total_runtime_seconds)
        .def("print", &Statistics::print)
    ;

}
