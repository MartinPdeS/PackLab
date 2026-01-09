#include "domain.h"

#include <pybind11/pybind11.h>
#include "pint/pint.h"

namespace py = pybind11;


PYBIND11_MODULE(interface_domain, module) {
    py::class_<Domain, std::shared_ptr<Domain>> domain_cls(module, "Domain", py::dynamic_attr());

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

                new (self.cast<Domain*>()) Domain(lx, ly, lz, use_periodic_boundaries);

                self.attr("_ureg") = ureg;  // Persist registry on the Python instance for later getters
            },
            py::arg("length_x"),
            py::arg("length_y"),
            py::arg("length_z"),
            py::arg("use_periodic_boundaries")
        )
        .def_property_readonly(
            "length_x",
            py::cpp_function(
                [](py::object self) {
                    const Domain& cpp_self = self.cast<const Domain&>();
                    py::object ureg = get_shared_ureg();
                    return meters_quantity_with_ureg(ureg, cpp_self.length_x);
                },
                py::is_method(domain_cls)
            ),
            "The length of the domain in the x dimension in meters."
        )
        .def_property_readonly(
            "length_y",
            py::cpp_function(
                [](py::object self) {
                    const Domain& cpp_self = self.cast<const Domain&>();
                    py::object ureg = get_shared_ureg();
                    return meters_quantity_with_ureg(ureg, cpp_self.length_y);
                },
                py::is_method(domain_cls)
            ),
            "The length of the domain in the y dimension in meters."
        )
        .def_property_readonly(
            "length_z",
            py::cpp_function(
                [](py::object self) {
                    const Domain& cpp_self = self.cast<const Domain&>();
                    py::object ureg = get_shared_ureg();
                    return meters_quantity_with_ureg(ureg, cpp_self.length_z);
                },
                py::is_method(domain_cls)
            ),
            "The length of the domain in the z dimension in meters."
        )
        .def_readonly(
            "use_periodic_boundaries",
            &Domain::use_periodic_boundaries,
            "Whether periodic boundary conditions are used in the domain."
        )
        .def_readonly(
            "volume",
            &Domain::volume,
            "The volume of the domain in cubic meters."
        )
        .def(
            "scale",
            &Domain::scale,
            py::arg("scale_factor"),
            "Scale the domain dimensions by the given factor."
        );
}
