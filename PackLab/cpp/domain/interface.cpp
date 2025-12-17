#include "domain.h"

#include <pybind11/pybind11.h>


PYBIND11_MODULE(interface_domain, module) {
    pybind11::class_<Domain>(module, "Domain")
        .def(
            pybind11::init<double, double, double, bool>(),
            pybind11::arg("length_x"),
            pybind11::arg("length_y"),
            pybind11::arg("length_z"),
            pybind11::arg("use_periodic_boundaries")
        )
        .def_readonly("length_x", &Domain::length_x)
        .def_readonly("length_y", &Domain::length_y)
        .def_readonly("length_z", &Domain::length_z)
        .def_readonly("use_periodic_boundaries", &Domain::use_periodic_boundaries)
        .def_readonly("volume", &Domain::volume)
        .def(
            "scale",
            &Domain::scale,
            pybind11::arg("scale_factor"),
            "Scale the domain dimensions by the given factor."
        )
    ;
}
