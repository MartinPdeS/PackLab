#include "result.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>

#include "pint/pint.h"

namespace py = pybind11;

PYBIND11_MODULE(result, module) {
    module.doc() = "Result class containing particle configuration and correlation functions";

    py::class_<Result, std::shared_ptr<Result>>(module, "PackingResult", R"doc(
Native result of a completed random sequential addition simulation.

Attributes
----------
sphere_configuration : PackingConfiguration
    Positions, radii, and class indices of all accepted spheres.
statistics : PackingStatistics
    Insertion counters, packing fractions, and runtime.
partial_volume_fractions : numpy.ndarray
    Volume fraction contributed by each radius class.

Notes
-----
The high-level :class:`PackLab.monte_carlo.PackingResult` adds plotting
and unit-aware convenience methods around this native object.
)doc")

        // -----------------------
        // Property: radii (N)
        // -----------------------
        .def_property_readonly(
            "radii",
            [](const Result& r) {const auto& radii = r.sphere_configuration->radii_values; return py::array_t<double>(radii.size(), radii.data());},
            "Particle radii as a NumPy array of shape (N)."
        )
        .def_readonly(
            "statistics",
            &Result::statistics,
            "Statistics about the simulation."
        )
        .def_readonly(
            "sphere_configuration",
            &Result::sphere_configuration,
            "Sphere configuration containing positions, radii, and classes."
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
            [](const Result& r) {const auto& classes = r.sphere_configuration->class_index_values; return py::array_t<int>(classes.size(), classes.data());},
            "Particle class index per particle (shape = (N))."
        )
        .def_readonly(
            "number_of_classes",
            &Result::number_of_classes,
            "Number of radius classes used to compute partial correlations."
        )
        // -----------------------
        // MCDomain
        // -----------------------
        .def_readonly(
            "domain",
            &Result::domain,
            py::return_value_policy::reference_internal,
            "Simulation domain object."
        )
        .def(
            "compute_partial_pair_distances",
            [](const Result& r, std::size_t maximum_pairs, std::uint64_t random_seed)
            {
                auto [distances, class_i, class_j] = r.compute_partial_pair_distances(maximum_pairs, random_seed);

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
            py::arg("random_seed") = 0,
            R"doc(
Sample unordered pair distances uniformly without replacement.

Parameters
----------
maximum_pairs : int
    Maximum number of pairs examined before radial filtering. Zero examines
    all pairs. Only distances below half the shortest box length are returned.
random_seed : int, default=0
    Unsigned seed for capped sampling; zero is deterministic.

Returns
-------
distances : numpy.ndarray
    Retained centre distances in SI metres, using minimum images when periodic.
class_i, class_j : numpy.ndarray
    Class labels of the two particles in each retained pair. The arrays may
    contain fewer than maximum_pairs entries because of radial filtering.
)doc"
        )
        // -----------------------
        // Partial g_ij(r)
        // -----------------------
        .def(
            "compute_partial_pair_correlation_function",
            [](const Result& r, std::size_t n_bins, std::size_t maximum_pairs, std::uint64_t random_seed)
            {
                auto [centers, gij] = r.compute_partial_pair_correlation_function(n_bins, maximum_pairs, random_seed);

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
            py::arg("random_seed") = 0,
            R"doc(
Estimate partial pair correlations from unordered particle pairs.

Parameters
----------
n_bins : int
    Positive number of radial bins up to half the shortest box length.
maximum_pairs : int, default=1000000
    Maximum number of pairs examined, sampled uniformly without replacement.
    Zero examines every pair. Pairs outside the radial range still count
    toward this limit.
random_seed : int, default=0
    Unsigned seed for capped pair sampling. Zero is deterministic.

Returns
-------
centers : numpy.ndarray
    Radial bin centres in SI metres.
g_ij : numpy.ndarray, shape (species, species, n_bins)
    Partial pair correlations with finite-particle normalization.
)doc"
        )
        .def("__repr__", [](const Result& self) {
            return "<PackingResult spheres=" + std::to_string(self.sphere_configuration->radii_values.size())
                + ", classes=" + std::to_string(self.number_of_classes) + ">";
        });
}
