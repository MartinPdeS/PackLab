#include "rsa.h"


#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>


PYBIND11_MODULE(interface_utils, module) {
    module.doc() = "Random sequential addition of non overlapping spheres in a 3D box";

    pybind11::class_<Vector3d>(module, "Vector3d")
        .def(pybind11::init<>())
        .def(pybind11::init<double, double, double>(), pybind11::arg("x"), pybind11::arg("y"), pybind11::arg("z"))
        .def_readwrite("x", &Vector3d::x)
        .def_readwrite("y", &Vector3d::y)
        .def_readwrite("z", &Vector3d::z)
    ;
}
