#pragma once

#include <mutex>
#include <stdexcept>
#include <vector>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class UnitRegistrySingleton {
public:
    static UnitRegistrySingleton& instance(); // declaration only

    UnitRegistrySingleton(const UnitRegistrySingleton&) = delete;
    UnitRegistrySingleton& operator=(const UnitRegistrySingleton&) = delete;

    void set_ureg(py::object ureg_object) {
        py::gil_scoped_acquire gil;

        std::call_once(initialization_flag, [&]() {
            ureg = std::move(ureg_object);
            if (ureg.is_none()) {
                throw std::runtime_error("UnitRegistrySingleton.set_ureg: ureg is None.");
            }
        });
    }

    py::object get_ureg() const {
        py::gil_scoped_acquire gil;

        if (ureg.is_none()) {
            throw std::runtime_error("UnitRegistrySingleton.get_ureg: ureg not initialized.");
        }
        return ureg;
    }

    bool is_initialized() const {
        py::gil_scoped_acquire gil;
        return !ureg.is_none();
    }

private:
    UnitRegistrySingleton() = default;

    mutable std::once_flag initialization_flag;
    py::object ureg = py::none();
};

py::object get_shared_ureg();

double to_meters_strict(py::handle value);

py::object registry_from_quantity(py::handle value);

py::object meters_quantity_with_ureg(const py::object& ureg, double meters_value);

std::vector<double> to_vector_units(py::handle values, const std::string& unit);

