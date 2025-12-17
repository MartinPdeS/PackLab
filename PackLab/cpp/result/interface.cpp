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
            [](const Result& r) {const auto& radii = r.particle_radii(); return py::array_t<double>(radii.size(), radii.data());},
            "Particle radii as a NumPy array of shape (N)."
        )
        .def_readonly(
            "statistics",
            &Result::statistics,
            "Statistics about the simulation."
        )
        .def_readonly(
            "partial_volume_fractions",
            &Result::partial_volume_fractions,
            "Partial volume fractions per class."
        )
        .def_readonly(
            "partial_volumes",
            &Result::partial_volumes,
            "Partial volumes per class."
        )
        // -----------------------
        // Property: particle classes (N)
        // -----------------------
        .def_property_readonly(
            "classes",
            [](const Result& r) {const auto& classes = r.sphere_configuration.class_index_values; return py::array_t<int>(classes.size(), classes.data());},
            "Particle class index per particle (shape = (N))."
        )
        .def_readonly(
            "number_of_classes",
            &Result::number_of_classes,
            "Number of radius classes used to compute partial correlations."
        )
        // -----------------------
        // Domain
        // -----------------------
        .def_readonly(
            "domain",
            &Result::domain,
            py::return_value_policy::reference_internal,
            "Simulation domain object."
        )
        .def(
            "compute_partial_pair_distances",
            [](const Result& r, std::size_t maximum_pairs)
            {
                auto [distances, class_i, class_j] = r.compute_partial_pair_distances(maximum_pairs);

                py::array::ShapeContainer shape({distances.size()});

                py::array_t<double> py_distance(shape);
                py::array_t<size_t> py_class_i(shape);
                py::array_t<size_t> py_class_j(shape);

                auto py_buffer_distance = py_distance.mutable_unchecked<1>();
                auto py_buffer_class_i = py_class_i.mutable_unchecked<1>();
                auto py_buffer_class_j = py_class_j.mutable_unchecked<1>();

                // Fill array
                for (std::size_t idx = 0; idx < distances.size(); ++idx)
                {
                    py_buffer_distance(idx) = distances[idx];
                    py_buffer_class_i(idx) = class_i[idx];
                    py_buffer_class_j(idx) = class_j[idx];
                }

                return py::make_tuple(py_distance, py_class_i, py_class_j);
            },
            py::arg("maximum_pairs"),
            "Compute all pair distances along with their class indices."
        )
        // -----------------------
        // Partial g_ij(r)
        // -----------------------
        .def(
            "compute_partial_pair_correlation_function",
            [](const Result& r, std::size_t n_bins, std::size_t maximum_pairs)
            {
                auto [centers, gij] = r.compute_partial_pair_correlation_function(n_bins, maximum_pairs);

                py::array::ShapeContainer shape({r.number_of_classes, r.number_of_classes, n_bins});

                // Allocate a NumPy array with shape (K, K, bins)
                py::array_t<double> gij_array(shape);

                // Access buffer
                auto buffer = gij_array.mutable_unchecked<3>();
                py::array_t<double> centers_array(centers.size(), centers.data());

                // Fill array
                for (std::size_t i = 0; i < gij.size(); ++i)
                    for (std::size_t j = 0; j < gij.size(); ++j)
                        for (std::size_t b = 0; b < n_bins; ++b)
                            buffer(i, j, b) = gij[i][j][b];


                return py::make_tuple(centers_array, gij_array);
            },
            py::arg("n_bins"),
            py::arg("maximum_pairs") = 1'000'000,
            "Compute partial pair correlation function g_ij(r) and return (centers, g_ij)."
        )
    ;
}
