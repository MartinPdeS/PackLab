#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <string>
#include <vector>

#include "analytical/domain_analytical/domain_analytical.h"
#include "pint/pint.h"
#include <iostream>

namespace py = pybind11;

static double quantity_scalar_to_meters(py::object quantity) {
    py::object ureg = get_shared_ureg();
    py::object meter = ureg.attr("meter");
    py::object converted = quantity.attr("to")(meter);
    return py::float_(converted.attr("magnitude"));
}

static std::vector<double> quantity_1d_to_meters_vector(py::object quantity) {
    py::object ureg = get_shared_ureg();
    py::object meter = ureg.attr("meter");

    py::object converted = quantity.attr("to")(meter);
    py::object magnitude = converted.attr("magnitude");

    py::module_ numpy = py::module_::import("numpy");
    py::object float64 = numpy.attr("float64");

    py::array magnitude_array = py::array::ensure(numpy.attr("asarray")(magnitude, float64));
    if (!magnitude_array || magnitude_array.ndim() != 1) {
        throw py::value_error("particle_radii must be a one dimensional quantity array.");
    }

    py::array_t<double, py::array::c_style | py::array::forcecast> radii_double(magnitude_array);
    auto buf = radii_double.unchecked<1>();

    std::vector<double> out(static_cast<std::size_t>(buf.shape(0)));
    for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
        out[static_cast<std::size_t>(i)] = buf(i);
    }
    return out;
}

static std::vector<double> array_like_1d_to_double_vector(py::object values) {
    py::module_ numpy = py::module_::import("numpy");
    py::object float64 = numpy.attr("float64");

    py::array array = py::array::ensure(numpy.attr("asarray")(values, float64));
    if (!array || array.ndim() != 1) {
        throw py::value_error("number_fractions must be a one dimensional array.");
    }

    py::array_t<double, py::array::c_style | py::array::forcecast> as_double(array);
    auto buf = as_double.unchecked<1>();

    std::vector<double> out(static_cast<std::size_t>(buf.shape(0)));
    for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
        out[static_cast<std::size_t>(i)] = buf(i);
    }
    return out;
}

static py::object meters_array_to_quantity(py::array_t<double> values_meters) {
    py::object ureg = get_shared_ureg();
    return values_meters * ureg.attr("meter");
}

static py::object meters3_array_to_quantity(py::array_t<double> values) {
    py::object ureg = get_shared_ureg();
    return values * (ureg.attr("meter**3"));
}

PYBIND11_MODULE(interface_domain_analytical, module) {
    module.doc() = "Polydisperse cubic domain with units handled in the pybind wrapper.";

    py::enum_<PolydisperseDomain::RoundingMode>(module, "RoundingMode")
        .value("floor", PolydisperseDomain::RoundingMode::Floor)
        .value("round", PolydisperseDomain::RoundingMode::Round);

    py::class_<PolydisperseDomain, std::shared_ptr<PolydisperseDomain>>(module, "Domain", py::dynamic_attr())
        .def(
            py::init([](
                py::object size_py,
                py::object particle_radii_py,
                py::object volume_fraction_py,
                py::object number_fractions_py,
                PolydisperseDomain::RoundingMode rounding_mode
            ) {
                const double size_m = quantity_scalar_to_meters(size_py);
                const std::vector<double> radii_m = quantity_1d_to_meters_vector(particle_radii_py);

                const double volume_fraction = py::float_(volume_fraction_py);  // dimensionless
                const std::vector<double> number_fraction = array_like_1d_to_double_vector(number_fractions_py);
                std::cout<<"number_fraction size: "<<number_fraction.size()<<std::endl;
                return PolydisperseDomain(size_m, radii_m, volume_fraction, number_fraction, rounding_mode);
            }),
            py::arg("size"),
            py::arg("particle_radii"),
            py::arg("volume_fraction"),
            py::arg("number_fractions"),
            py::arg("rounding_mode") = PolydisperseDomain::RoundingMode::Floor
        )
        .def_property_readonly(
            "size",
            [](const PolydisperseDomain& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.size_meters) * ureg.attr("meter");
            }
        )
        .def_property_readonly(
            "particle_radii",
            [](const PolydisperseDomain& self) {
                py::array_t<double> arr(static_cast<py::ssize_t>(self.particle_radii.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = self.particle_radii[static_cast<std::size_t>(i)];
                }
                return meters_array_to_quantity(arr);
            }
        )
        .def_readonly(
            "volume_fraction",
            &PolydisperseDomain::volume_fraction
        )
        .def_property_readonly(
            "number_fractions",
            [](const PolydisperseDomain& self) {
                py::array_t<double> arr(static_cast<py::ssize_t>(self.number_fractions.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = self.number_fractions[static_cast<std::size_t>(i)];
                }
                return arr;
            }
        )
        .def_property_readonly(
            "volume",
            [](const PolydisperseDomain& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.volume()) * (ureg.attr("meter**3"));
            }
        )
        .def_property_readonly(
            "particle_volumes",
            [](const PolydisperseDomain& self) {
                const auto values = self.get_particle_volumes();
                py::array_t<double> arr(static_cast<py::ssize_t>(values.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = values[static_cast<std::size_t>(i)];
                }
                return meters3_array_to_quantity(arr);
            }
        )
        .def_property_readonly(
            "total_particle_volume",
            [](const PolydisperseDomain& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_total_particle_volume()) * (ureg.attr("meter**3"));
            }
        )
        .def_property_readonly(
            "mean_particle_volume_number_weighted",
            [](const PolydisperseDomain& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_mean_particle_volume_number_weighted()) * (ureg.attr("meter**3"));
            }
        )
        .def_property_readonly(
            "number_of_particles_total",
            &PolydisperseDomain::get_number_of_particles_total
        )
        .def_property_readonly(
            "number_of_particles_per_radius",
            [](const PolydisperseDomain& self) {
                const auto values = self.get_number_of_particles_per_radius();
                py::array_t<long long> arr(static_cast<py::ssize_t>(values.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = static_cast<long long>(values[static_cast<std::size_t>(i)]);
                }
                return arr;
            }
        )
        .def_property_readonly(
            "particle_densities_per_radius",
            [](const PolydisperseDomain& self) {
                const auto values = self.get_particle_densities_per_radius();
                py::array_t<double> arr(static_cast<py::ssize_t>(values.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = values[static_cast<std::size_t>(i)];
                }
                py::object ureg = get_shared_ureg();
                return arr * (py::float_(1.0) / (ureg.attr("meter**3")));
            }
        )
        .def_property_readonly(
            "particle_density_total",
            [](const PolydisperseDomain& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_particle_density_total()) * (py::float_(1.0) / (ureg.attr("meter**3")));
            }
        )
        .def_property_readonly(
            "volume_fraction_per_radius",
            [](const PolydisperseDomain& self) {
                const auto values = self.get_volume_fraction_per_radius();
                py::array_t<double> arr(static_cast<py::ssize_t>(values.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = values[static_cast<std::size_t>(i)];
                }
                return arr;
            }
        )
        .def(
            "sample_particle_radii",
            [](const PolydisperseDomain& self, long long number_of_samples, py::object seed_py) {
                std::uint64_t seed = 0;
                if (seed_py.is_none())
                    seed = static_cast<std::uint64_t>(std::random_device{}());
                else
                    seed = py::int_(seed_py);


                const auto samples = self.sample_particle_radii(number_of_samples, seed);

                py::array_t<double> arr(static_cast<py::ssize_t>(samples.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = samples[static_cast<std::size_t>(i)];
                }
                return meters_array_to_quantity(arr);
            },
            py::arg("number_of_samples"),
            py::arg("seed") = py::none()
        )
        .def(
            "print_bins",
            // &PolydisperseDomain::print_bins,
            [](const PolydisperseDomain& self, int precision){self.print_bins(precision);},
            py::arg("precision") = 6,
            "Print a per bin table in base units."
        )
        .def(
            "bins_table",
            &PolydisperseDomain::bins_table,
            py::arg("precision") = 6,
            "Return the per bin table as a string in base units."
        )

        ;
}
