#include "statistics.h"

#include <pybind11/pybind11.h>
#include <sstream>


PYBIND11_MODULE(statistics, module) {
    pybind11::class_<Statistics>(module, "PackingStatistics", R"doc(
Summary statistics from an RSA simulation.

Attributes
----------
sphere_count : int
    Number of accepted spheres.
packing_fraction_geometry : float
    Packing fraction calculated from accepted sphere volumes.
total_runtime_seconds : float
    Wall-clock runtime of the simulation.
)doc")
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
        .def("print", &Statistics::print, "Print a human-readable statistics summary.")
        .def("__repr__", [](const Statistics& self) {
            std::ostringstream stream;
            stream << "<PackingStatistics spheres=" << self.sphere_count
                   << ", packing_fraction=" << self.packing_fraction_geometry << ">";
            return stream.str();
        });

}
