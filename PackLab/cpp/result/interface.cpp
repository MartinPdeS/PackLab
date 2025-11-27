#include "result.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


PYBIND11_MODULE(interface_result, module) {
    pybind11::class_<Result>(module, "Result")
        .def(
            pybind11::init<std::vector<Vector3d>, std::vector<double>, Domain>(),
            pybind11::arg("positions"),
            pybind11::arg("radii"),
            pybind11::arg("domain")
        )
        .def(
            pybind11::init<std::vector<std::vector<double>>, std::vector<double>, Domain>(),
            pybind11::arg("positions"),
            pybind11::arg("radii"),
            pybind11::arg("domain")
        )
        .def_property_readonly(
            "positions",
            [](const Result& r) {
                const auto& positions = r.particle_positions();
                pybind11::array_t<double> array({positions.size(), static_cast<std::size_t>(3)});
                auto buffer = array.mutable_unchecked<2>();
                for (std::size_t i = 0; i < positions.size(); ++i) {
                    buffer(i, 0) = positions[i].x;
                    buffer(i, 1) = positions[i].y;
                    buffer(i, 2) = positions[i].z;
                }
                return array;
            },
            "Particle positions as a NumPy array of shape (N, 3)"
        )

        .def_property_readonly(
            "radii",
            [](const Result& r) {
                const auto& radii = r.particle_radii();
                return pybind11::array_t<double>(radii.size(), radii.data());
            },
            "Particle radii as a NumPy array of shape (N)"
        )

        .def_property_readonly(
            "domain_box",
            [](const Result& r) -> const Domain& {
                return r.domain;
            },
            pybind11::return_value_policy::reference_internal,
            "Simulation domain box object"
        )

        .def_property_readonly(
            "pair_correlation_centers",
            [](const Result& r) {
                const auto& centers = r.pair_correlation_centers();
                return pybind11::array_t<double>(centers.size(), centers.data());
            },
            "Bin centers for the pair correlation function"
        )

        .def_property_readonly(
            "pair_correlation_values",
            [](const Result& r) {
                const auto& values = r.pair_correlation_values();
                return pybind11::array_t<double>(values.size(), values.data());
            },
            "Pair correlation values g(r)"
        )
        .def_readonly(
            "domain",
            &Result::domain,
            "Simulation domain box object"
        )

        // Methods
        .def(
            "compute_pair_correlation_function",
            &Result::compute_pair_correlation_function,
            pybind11::arg("bins") = 90,
            pybind11::arg("maximum_number_of_pairs") = 250000,
            pybind11::arg("random_seed") = 0,
            "Compute the pair correlation function and store the result internally"
        );

}
