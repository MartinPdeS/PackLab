#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <sstream>

#include "py_solver.h"
#include "pint/pint.h"

namespace py = pybind11;

static double quantity_scalar_to_double_in_units(py::object quantity, const char* unit_string) {
    py::object converted = quantity.attr("to")(py::str(unit_string));
    return py::float_(converted.attr("magnitude"));
}

static std::vector<double> quantity_1d_to_double_vector_in_units(py::object quantity, const char* unit_string) {
    py::object converted = quantity.attr("to")(py::str(unit_string));
    py::object magnitude = converted.attr("magnitude");

    py::module_ numpy = py::module_::import("numpy");
    py::object float64 = numpy.attr("float64");
    py::array arr = py::array::ensure(numpy.attr("asarray")(magnitude, float64));
    if (!arr || arr.ndim() != 1) {
        throw py::value_error("Expected a one dimensional quantity array.");
    }

    py::array_t<double, py::array::c_style | py::array::forcecast> arr_double(arr);
    auto buf = arr_double.unchecked<1>();

    std::vector<double> out(static_cast<std::size_t>(buf.shape(0)));
    for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
        out[static_cast<std::size_t>(i)] = buf(i);
    }
    return out;
}

static py::array_t<double> vector_to_numpy_1d(const std::vector<double>& v) {
    py::array_t<double> arr(static_cast<py::ssize_t>(v.size()));
    auto buf = arr.mutable_unchecked<1>();
    for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
        buf(i) = v[static_cast<std::size_t>(i)];
    }
    return arr;
}

static py::array_t<double> vector_to_numpy_2d_row_major(const std::vector<double>& v, std::size_t n0, std::size_t n1) {
    py::array_t<double> arr({static_cast<py::ssize_t>(n0), static_cast<py::ssize_t>(n1)});
    auto buf = arr.mutable_unchecked<2>();
    for (std::size_t i = 0; i < n0; ++i) {
        for (std::size_t j = 0; j < n1; ++j) {
            buf(static_cast<py::ssize_t>(i), static_cast<py::ssize_t>(j)) = v[i * n1 + j];
        }
    }
    return arr;
}

static py::array_t<double> vector_to_numpy_3d_row_major(const std::vector<double>& v, std::size_t n0, std::size_t n1, std::size_t n2) {
    py::array_t<double> arr({static_cast<py::ssize_t>(n0), static_cast<py::ssize_t>(n1), static_cast<py::ssize_t>(n2)});
    auto buf = arr.mutable_unchecked<3>();
    std::size_t index = 0;
    for (std::size_t i = 0; i < n0; ++i) {
        for (std::size_t j = 0; j < n1; ++j) {
            for (std::size_t k = 0; k < n2; ++k) {
                buf(static_cast<py::ssize_t>(i), static_cast<py::ssize_t>(j), static_cast<py::ssize_t>(k)) = v[index++];
            }
        }
    }
    return arr;
}

PYBIND11_MODULE(solver, module) {
    module.doc() = "Percus Yevick mixture solver (C++ core, Pint handled in wrapper).";

    py::class_<PercusYevickResult>(module, "PercusYevickResult", py::dynamic_attr(), R"doc(
Result of a multicomponent Percus-Yevick calculation.

Attributes
----------
radii : pint.Quantity, shape (n_species,)
    Particle radii in meters.
densities : pint.Quantity, shape (n_species,)
    Number densities in inverse cubic meters.
wavenumber : pint.Quantity, shape (n_wavenumber,)
    Wavenumber grid in inverse meters used for the inverse transform.
g, h, H : numpy.ndarray
    Real-space and reciprocal-space correlation functions.
)doc")
        .def_property_readonly("epsilons", [](const PercusYevickResult& r) {
            return vector_to_numpy_1d(r.epsilons);
        })
        .def_property_readonly("radii", [](const PercusYevickResult& r) {
            py::object ureg = get_shared_ureg();
            return vector_to_numpy_1d(r.radii_m) * ureg.attr("meter");
        })
        .def_property_readonly("densities", [](const PercusYevickResult& r) {
            py::object ureg = get_shared_ureg();
            py::object inv_m3 = py::float_(1.0) / ureg.attr("meter**3");
            return vector_to_numpy_1d(r.densities_per_m3) * inv_m3;
        })
        .def_property_readonly("wavenumber", [](const PercusYevickResult& r) {
            py::object ureg = get_shared_ureg();
            py::object inv_m = py::float_(1.0) / ureg.attr("meter");
            return vector_to_numpy_1d(r.wavenumber_per_m) * inv_m;
        })
        .def_property_readonly("distances", [](const PercusYevickResult& r) {
            py::object ureg = get_shared_ureg();
            return vector_to_numpy_1d(r.distances_m) * ureg.attr("meter");
        })
        .def_property_readonly("R_ij", [](const PercusYevickResult& r) {
            py::object ureg = get_shared_ureg();
            return vector_to_numpy_2d_row_major(r.R_ij_m, r.number_of_species, r.number_of_species) * ureg.attr("meter");
        })
        .def_property_readonly("S_ij", [](const PercusYevickResult& r) {
            py::object ureg = get_shared_ureg();
            return vector_to_numpy_2d_row_major(r.S_ij_m, r.number_of_species, r.number_of_species) * ureg.attr("meter");
        })
        .def_property_readonly("A_i", [](const PercusYevickResult& r) {
            return vector_to_numpy_1d(r.A_i);
        })
        .def_property_readonly("B_i", [](const PercusYevickResult& r) {
            py::object ureg = get_shared_ureg();
            return vector_to_numpy_1d(r.B_i_m) * ureg.attr("meter");
        })
        .def_property_readonly("D_ij", [](const PercusYevickResult& r) {
            py::object ureg = get_shared_ureg();
            return vector_to_numpy_2d_row_major(r.D_ij_m2, r.number_of_species, r.number_of_species) * ureg.attr("meter**2");
        })
        .def_property_readonly("Cpy", [](const PercusYevickResult& r) {
            return vector_to_numpy_3d_row_major(r.Cpy, r.number_of_species, r.number_of_species, r.number_of_wavenumber_points);
        })
        .def_property_readonly("H", [](const PercusYevickResult& r) {
            return vector_to_numpy_3d_row_major(r.H, r.number_of_species, r.number_of_species, r.number_of_wavenumber_points);
        })
        .def_property_readonly("h", [](const PercusYevickResult& r) {
            return vector_to_numpy_3d_row_major(r.h, r.number_of_species, r.number_of_species, r.number_of_r_points);
        })
        .def_property_readonly("g", [](const PercusYevickResult& r) {
            return vector_to_numpy_3d_row_major(r.g, r.number_of_species, r.number_of_species, r.number_of_r_points);
        })
        .def("__repr__", [](const PercusYevickResult& self) {
            return "<PercusYevickResult species=" + std::to_string(self.number_of_species)
                + ", wavenumber_points=" + std::to_string(self.number_of_wavenumber_points) + ">";
        });

    py::class_<PercusYevickSolver, std::shared_ptr<PercusYevickSolver>>(module, "PercusYevickSolver", R"doc(
Solve the multicomponent Percus-Yevick hard-sphere model.

Parameters
----------
densities : pint.Quantity, shape (n_species,)
    Species number densities, convertible to inverse cubic meters.
radii : pint.Quantity, shape (n_species,)
    Species radii, convertible to meters.
wavenumber : pint.Quantity or "auto", default="auto"
    Strictly increasing wavenumber grid, convertible to inverse meters. ``"auto"``
    derives a grid when :meth:`compute` receives the requested distances.
radial_resolution : pint.Quantity or None, optional
    Smallest real-space feature to resolve when ``wavenumber="auto"``. By
    default this is one twentieth of the smallest particle radius.
samples_per_oscillation : int, default=12
    Reciprocal-space samples per sinc-kernel oscillation at the largest
    requested distance when ``wavenumber="auto"``. Must be at least 8.

Notes
-----
The inverse radial transform must resolve the sinc kernel at the largest
requested distance. :meth:`compute` emits a ``RuntimeWarning`` when this grid
has fewer than eight samples per kernel oscillation. Use
``PackLab.analytical.make_wavenumber_grid`` to construct a balanced grid.
)doc")
        .def(
            py::init([](
                py::object densities_py,
                py::object radii_py,
                py::object wavenumber_py,
                py::object radial_resolution_py,
                int samples_per_oscillation
            ) {
                const std::vector<double> densities_per_m3 = quantity_1d_to_double_vector_in_units(densities_py, "1/meter**3");
                const std::vector<double> radii_m = quantity_1d_to_double_vector_in_units(radii_py, "meter");

                if (py::isinstance<py::str>(wavenumber_py)) {
                    const std::string mode = py::cast<std::string>(wavenumber_py);
                    if (mode != "auto") {
                        throw py::value_error("wavenumber must be a quantity array or 'auto'.");
                    }
                    const double radial_resolution_m = radial_resolution_py.is_none()
                        ? 0.0
                        : quantity_scalar_to_double_in_units(radial_resolution_py, "meter");
                    return std::make_shared<PercusYevickSolver>(
                        densities_per_m3,
                        radii_m,
                        radial_resolution_m,
                        samples_per_oscillation
                    );
                }

                if (!radial_resolution_py.is_none() || samples_per_oscillation != 12) {
                    throw py::value_error(
                        "radial_resolution and samples_per_oscillation are only available when wavenumber='auto'."
                    );
                }
                const std::vector<double> wavenumber_per_m = quantity_1d_to_double_vector_in_units(wavenumber_py, "1/meter");
                return std::make_shared<PercusYevickSolver>(densities_per_m3, radii_m, wavenumber_per_m);
            }),
            py::arg("densities"),
            py::arg("radii"),
            py::arg("wavenumber") = py::str("auto"),
            py::arg("radial_resolution") = py::none(),
            py::arg("samples_per_oscillation") = 12,
            "Create a Percus-Yevick solver."
        )
        .def(
            "compute",
            [](const PercusYevickSolver& self, py::object distances_py) {
                const std::vector<double> distances_m = quantity_1d_to_double_vector_in_units(distances_py, "meter");
                if (!self.uses_automatic_wavenumber()) {
                    if (const auto warning = self.radial_grid_warning(distances_m)) {
                        if (PyErr_WarnEx(PyExc_RuntimeWarning, warning->c_str(), 2) != 0) {
                            throw py::error_already_set();
                        }
                    }
                }
                return self.compute(distances_m);
            },
            py::arg("distances"),
            R"doc(
Compute the correlation functions on a real-space distance grid.

Parameters
----------
distances : pint.Quantity, shape (n_r,)
    Radial positions, convertible to meters.

Returns
-------
PercusYevickResult
    Computed structural and correlation functions.
)doc"
        )
        .def("__repr__", [](const PercusYevickSolver& self) {
            return self.uses_automatic_wavenumber()
                ? "<PercusYevickSolver wavenumber='auto'>"
                : "<PercusYevickSolver>";
        });
}
