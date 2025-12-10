#include "result.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(interface_result, module) {
    module.doc() = "Result class containing particle configuration and correlation functions";

    py::class_<Result>(module, "Result")

        // -----------------------
        // Property: positions (N,3)
        // -----------------------
        .def_property_readonly(
            "positions",
            [](const Result& r) {
                const auto& positions = r.particle_positions();
                py::array_t<double> array({positions.size(), (std::size_t)3});
                auto buf = array.mutable_unchecked<2>();
                for (std::size_t i = 0; i < positions.size(); ++i) {
                    buf(i, 0) = positions[i].x;
                    buf(i, 1) = positions[i].y;
                    buf(i, 2) = positions[i].z;
                }
                return array;
            },
            "Particle positions as a NumPy array of shape (N, 3)."
        )

        // -----------------------
        // Property: radii (N)
        // -----------------------
        .def_property_readonly(
            "radii",
            [](const Result& r) {
                const auto& radii = r.particle_radii();
                return py::array_t<double>(radii.size(), radii.data());
            },
            "Particle radii as a NumPy array of shape (N)."
        )

        .def_readonly(
            "statistics",
            &Result::statistics,
            "Statistics about the simulation."
        )

        // -----------------------
        // Property: particle classes (N)
        // -----------------------
        .def_property_readonly(
            "classes",
            [](const Result& r) {
                const auto& classes = r.sphere_configuration.class_index_values_;
                return py::array_t<int>(classes.size(), classes.data());
            },
            "Particle class index per particle (shape = (N))."
        )

        .def_readonly(
            "number_of_classes",
            &Result::number_of_classes_,
            "Number of radius classes used to compute partial correlations."
        )

        // -----------------------
        // Domain
        // -----------------------
        .def_property_readonly(
            "domain",
            [](const Result& r) -> const Domain& {
                return r.domain;
            },
            py::return_value_policy::reference_internal,
            "Simulation domain object."
        )

        // -----------------------
        // Pair correlation outputs
        // -----------------------
        .def_property_readonly(
            "pair_correlation_centers",
            [](const Result& r) {
                const auto& centers = r.pair_correlation_centers();
                return py::array_t<double>(centers.size(), centers.data());
            },
            "g(r) radial bin centers."
        )
        .def_property_readonly(
            "pair_correlation_values",
            [](const Result& r) {
                const auto& values = r.pair_correlation_values();
                return py::array_t<double>(values.size(), values.data());
            },
            "Computed g(r) values."
        )
        .def_property_readonly(
            "pair_correlation_mean_values",
            [](const Result& r) {
                const auto& v = r.pair_correlation_mean_values();
                return py::array_t<double>(v.size(), v.data());
            },
            "Mean g(r) values (when multiple trials are aggregated)."
        )
        .def_property_readonly(
            "pair_correlation_std_values",
            [](const Result& r) {
                const auto& v = r.pair_correlation_std_values();
                return py::array_t<double>(v.size(), v.data());
            },
            "Standard deviation of g(r) values."
        )

        // -----------------------
        // Compute g(r)
        // -----------------------
        .def(
            "compute_pair_correlation_function",
            &Result::compute_pair_correlation_function,
            py::arg("bins") = 90,
            py::arg("maximum_distance") = 0.0,
            py::arg("method") = "grid",
            py::arg("maximum_number_of_pairs") = 100000000,
            py::arg("random_seed") = 0,
            py::arg("grid_cell_size") = 0.0,
            "Compute the pair correlation function and store results inside this object."
        )

        // -----------------------
        // Partial g_ij(r)
        // -----------------------
        .def(
            "compute_partial_pair_correlation_function",
            [](const Result& r,
            std::size_t bins,
            std::size_t maximum_pairs,
            std::uint64_t random_seed)
            {
                auto [centers, gij] =
                    r.compute_partial_pair_correlation_function(
                        bins,
                        maximum_pairs,
                        random_seed
                    );

                const int K = static_cast<int>(gij.size());

                // Allocate a NumPy array with shape (K, K, bins)
                py::array_t<double> gij_array(
                    py::array::ShapeContainer{
                        static_cast<py::ssize_t>(K),
                        static_cast<py::ssize_t>(K),
                        static_cast<py::ssize_t>(bins)
                    }
                );

                // Access buffer
                auto buf = gij_array.mutable_unchecked<3>();

                // Fill array
                for (int i = 0; i < K; ++i)
                    for (int j = 0; j < K; ++j)
                        for (std::size_t b = 0; b < bins; ++b)
                            buf(i, j, static_cast<py::ssize_t>(b)) = gij[i][j][b];

                // Return (centers, gij_array)
                py::array_t<double> centers_array(
                    centers.size(),
                    centers.data()
                );

                return py::make_tuple(
                    centers_array,
                    gij_array
                );
            },
            py::arg("number_of_distance_bins"),
            py::arg("maximum_pairs") = 1'000'000,
            py::arg("random_seed") = 0,
            "Compute partial pair correlation function g_ij(r) and return (centers, g_ij)."
        )
    ;
}
