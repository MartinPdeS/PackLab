#include "domain.h"

#include <pybind11/pybind11.h>
#include <sstream>
#include "pint/pint.h"

namespace py = pybind11;


PYBIND11_MODULE(domain, module) {
    py::class_<MCDomain, std::shared_ptr<MCDomain>> domain_cls(module, "PackingDomain", py::dynamic_attr(), R"doc(
Three-dimensional rectangular domain for random sequential addition.

Parameters
----------
length_x, length_y, length_z : pint.Quantity
    Positive box lengths, convertible to meters.
use_periodic_boundaries : bool
    Whether to apply periodic boundary conditions on all three axes.

Attributes
----------
length_x, length_y, length_z : pint.Quantity
    Box lengths in meters.
volume : float
    Box volume in cubic meters.
)doc");

    domain_cls
        .def(
            "__init__",
            [](
                py::object self,
                py::object length_x,
                py::object length_y,
                py::object length_z,
                bool use_periodic_boundaries
            ) {
                py::object ureg = get_shared_ureg();

                const double lx = to_meters_strict(length_x);
                const double ly = to_meters_strict(length_y);
                const double lz = to_meters_strict(length_z);

                new (self.cast<MCDomain*>()) MCDomain(lx, ly, lz, use_periodic_boundaries);

                self.attr("_ureg") = ureg;  // Persist registry on the Python instance for later getters
            },
            py::arg("length_x"),
            py::arg("length_y"),
            py::arg("length_z"),
            py::arg("use_periodic_boundaries"),
            "Create a packing domain from unit-bearing box lengths."
        )
        .def_property_readonly(
            "length_x",
            py::cpp_function(
                [](py::object self) {
                    const MCDomain& cpp_self = self.cast<const MCDomain&>();
                    py::object ureg = get_shared_ureg();
                    return meters_quantity_with_ureg(ureg, cpp_self.length_x);
                },
                py::is_method(domain_cls)
            ),
            "pint.Quantity: Length of the x axis in meters."
        )
        .def_property_readonly(
            "length_y",
            py::cpp_function(
                [](py::object self) {
                    const MCDomain& cpp_self = self.cast<const MCDomain&>();
                    py::object ureg = get_shared_ureg();
                    return meters_quantity_with_ureg(ureg, cpp_self.length_y);
                },
                py::is_method(domain_cls)
            ),
            "pint.Quantity: Length of the y axis in meters."
        )
        .def_property_readonly(
            "length_z",
            py::cpp_function(
                [](py::object self) {
                    const MCDomain& cpp_self = self.cast<const MCDomain&>();
                    py::object ureg = get_shared_ureg();
                    return meters_quantity_with_ureg(ureg, cpp_self.length_z);
                },
                py::is_method(domain_cls)
            ),
            "pint.Quantity: Length of the z axis in meters."
        )
        .def_readonly(
            "use_periodic_boundaries",
            &MCDomain::use_periodic_boundaries,
            "bool: Whether periodic boundary conditions are enabled."
        )
        .def_readonly(
            "volume",
            &MCDomain::volume,
            "float: Box volume in cubic meters."
        )
        .def(
            "scale",
            &MCDomain::scale,
            py::arg("scale_factor"),
            R"doc(
Scale all box lengths in place.

Parameters
----------
scale_factor : float
    Positive multiplicative scale factor.
)doc"
        )
        .def("__repr__", [](const MCDomain& self) {
            std::ostringstream stream;
            stream << "<PackingDomain lengths=(" << self.length_x << ", " << self.length_y << ", "
                   << self.length_z << ") m, periodic=" << (self.use_periodic_boundaries ? "True" : "False") << ">";
            return stream.str();
        });
}
