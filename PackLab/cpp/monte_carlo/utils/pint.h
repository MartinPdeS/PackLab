#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>

namespace py = pybind11;

double to_meters_strict(py::handle value) {
    // No bare floats or ints allowed
    if (py::isinstance<py::float_>(value) || py::isinstance<py::int_>(value)) {
        throw py::type_error("length must be a pint.Quantity (numbers without units are not accepted)");
    }

    if (!py::hasattr(value, "to")) {
        throw py::type_error("length must be a pint.Quantity (must provide a .to method)");
    }

    py::object q = py::reinterpret_borrow<py::object>(value);
    py::object q_m = q.attr("to")("meter");

    if (!py::hasattr(q_m, "magnitude")) {
        throw py::type_error("length.to('meter') did not return an object with a magnitude attribute");
    }

    return py::cast<double>(q_m.attr("magnitude"));
}

py::object registry_from_quantity(py::handle value) {
    // Pint Quantity objects expose a registry as _REGISTRY
    if (!py::hasattr(value, "_REGISTRY")) {
        throw py::type_error("length must be a pint.Quantity (missing _REGISTRY)");
    }
    return py::reinterpret_borrow<py::object>(value).attr("_REGISTRY");
}

py::object meters_quantity_with_ureg(const py::object& ureg, double meters_value) {
    return ureg.attr("Quantity")(meters_value, "meter");
}


std::vector<double> to_vector_units(py::handle values, const std::string& unit) {
    if (!py::hasattr(values, "__iter__")) {
        throw py::type_error("expected an iterable of pint.Quantity objects");
    }

    std::vector<double> magnitudes;

    // Reserve only when __len__ exists to avoid surprising errors for generators
    if (py::hasattr(values, "__len__")) {
        magnitudes.reserve(py::len(values));
    }

    for (py::handle item : py::reinterpret_borrow<py::iterable>(values)) {
        if (py::isinstance<py::int_>(item) || py::isinstance<py::float_>(item)) {
            throw py::type_error("expected pint.Quantity elements (numbers without units are not accepted)");
        }
        if (!py::hasattr(item, "to")) {
            throw py::type_error("expected pint.Quantity elements (missing .to method)");
        }

        py::object quantity = py::reinterpret_borrow<py::object>(item);
        py::object converted = quantity.attr("to")(unit);

        if (!py::hasattr(converted, "magnitude")) {
            throw py::type_error("quantity.to(unit) did not return an object with a magnitude attribute");
        }

        magnitudes.push_back(py::cast<double>(converted.attr("magnitude")));
    }

    return magnitudes;
}
