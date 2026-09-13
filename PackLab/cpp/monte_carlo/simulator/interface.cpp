#include "simulator.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <cmath>
#include <sstream>

#include <pint/pint.h>

namespace py = pybind11;

static py::array_t<double> vector3d_list_to_numpy(const std::vector<Vector3d>& values) {
    const py::ssize_t count = static_cast<py::ssize_t>(values.size());
    auto array = py::array_t<double>({count, py::ssize_t(3)});
    auto buffer = array.mutable_unchecked<2>();

    for (py::ssize_t index = 0; index < count; ++index) {
        buffer(index, 0) = values[static_cast<std::size_t>(index)].x;
        buffer(index, 1) = values[static_cast<std::size_t>(index)].y;
        buffer(index, 2) = values[static_cast<std::size_t>(index)].z;
    }

    return array;
}

static py::array_t<double> double_list_to_numpy(const std::vector<double>& values) {
    const py::ssize_t count = static_cast<py::ssize_t>(values.size());
    auto array = py::array_t<double>(count);
    std::memcpy(array.mutable_data(), values.data(), static_cast<std::size_t>(count) * sizeof(double));
    return array;
}

static std::shared_ptr<SphereConfiguration> configuration_from_arrays(
    const py::object& positions_quantity,
    const py::object& radii_quantity,
    const py::object& classes
) {
    py::array class_array = py::array::ensure(classes);
    if (!class_array) {
        throw py::type_error("classes_index must be an integer array.");
    }
    py::module_ numpy = py::module_::import("numpy");
    if (!py::cast<bool>(numpy.attr("issubdtype")(class_array.attr("dtype"), numpy.attr("integer")))) {
        throw py::type_error("classes_index must be an integer array.");
    }
    py::array_t<double, py::array::c_style | py::array::forcecast> positions(
        positions_quantity.attr("to")("meter").attr("magnitude")
    );
    py::array_t<double, py::array::c_style | py::array::forcecast> radii(
        radii_quantity.attr("to")("meter").attr("magnitude")
    );
    py::array_t<int, py::array::c_style | py::array::forcecast> class_indices(class_array);

    if (positions.ndim() != 2 || positions.shape(1) != 3) {
        throw py::value_error("positions must have shape (N, 3).");
    }
    if (radii.ndim() != 1 || class_indices.ndim() != 1 ||
        radii.shape(0) != positions.shape(0) || class_indices.shape(0) != positions.shape(0)) {
        throw py::value_error("radii and classes_index must be one-dimensional with one value per position.");
    }

    auto configuration = std::make_shared<SphereConfiguration>();
    const auto position_values = positions.unchecked<2>();
    const auto radius_values = radii.unchecked<1>();
    const auto class_values = class_indices.unchecked<1>();
    for (py::ssize_t index = 0; index < positions.shape(0); ++index) {
        const double x = position_values(index, 0);
        const double y = position_values(index, 1);
        const double z = position_values(index, 2);
        const double radius = radius_values(index);
        const int class_index = class_values(index);
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
            throw py::value_error("positions must contain only finite values.");
        }
        if (!std::isfinite(radius) || radius <= 0.0) {
            throw py::value_error("radii must be finite and positive.");
        }
        if (class_index < 0) {
            throw py::value_error("classes_index must contain non-negative integers.");
        }
        configuration->center_positions.push_back({x, y, z});
        configuration->radii_values.push_back(radius);
        configuration->class_index_values.push_back(class_index);
    }
    return configuration;
}

py::object run_and_wrap(Simulator& self) {
    auto cpp_result = self.run();  // whatever your C++ returns

    py::object ResultClass = py::module_::import("PackLab.monte_carlo.results").attr("PackingResult");

    py::dict metadata;
    metadata["random_seed"] = self.get_options()->random_seed;
    metadata["maximum_attempts"] = self.get_options()->maximum_attempts;
    metadata["maximum_spheres"] = self.get_options()->maximum_spheres;
    metadata["target_packing_fraction"] = self.get_options()->target_packing_fraction;
    return ResultClass(
        py::arg("binding") = py::cast(std::move(cpp_result)),
        py::arg("source") = "rsa",
        py::arg("run_metadata") = metadata
    );
}

PYBIND11_MODULE(simulator, module) {
    module.doc() = "Random sequential addition of non overlapping spheres in a 3D box";

    py::class_<SphereConfiguration, std::shared_ptr<SphereConfiguration>>(module, "PackingConfiguration", R"doc(
Sphere centers, radii, and size-class labels of a packing.

)doc")
        .def_property_readonly(
            "count",
            [](const std::shared_ptr<SphereConfiguration> sphere_configuration) { return sphere_configuration->radii().size(); },
            "Number of spheres in the configuration."
        )
        .def_static(
            "from_arrays",
            &configuration_from_arrays,
            py::arg("positions"),
            py::arg("radii"),
            py::arg("classes_index"),
            R"doc(
Create a configuration from validated, unit-bearing particle arrays.

Parameters
----------
positions : pint.Quantity, shape (N, 3)
    Sphere centres convertible to metres.
radii : pint.Quantity, shape (N,)
    Positive sphere radii convertible to metres.
classes_index : array-like, shape (N,)
    Non-negative integer class label for each sphere.
)doc"
        )
        .def(
            "compute_partial_pair_correlation_function",
            [](const std::shared_ptr<SphereConfiguration>& configuration,
               const std::shared_ptr<MCDomain>& domain,
               std::size_t n_bins,
               std::size_t maximum_pairs) {
                if (configuration->center_positions.size() < 2) {
                    throw py::value_error("at least two particles are required for partial pair correlation.");
                }
                int maximum_class = -1;
                for (const int class_index : configuration->class_index_values) {
                    maximum_class = std::max(maximum_class, class_index);
                }
                Result result(configuration, domain, Statistics{}, static_cast<std::size_t>(maximum_class + 1));
                auto [centers, gij] = result.compute_partial_pair_correlation_function(n_bins, maximum_pairs);
                py::array_t<double> centers_array(centers.size(), centers.data());
                py::array_t<double> gij_array({
                    static_cast<py::ssize_t>(maximum_class + 1),
                    static_cast<py::ssize_t>(maximum_class + 1),
                    static_cast<py::ssize_t>(n_bins)
                });
                auto output = gij_array.mutable_unchecked<3>();
                for (std::size_t i = 0; i < gij.size(); ++i)
                    for (std::size_t j = 0; j < gij.size(); ++j)
                        for (std::size_t bin = 0; bin < n_bins; ++bin)
                            output(i, j, bin) = gij[i][j][bin];
                return py::make_tuple(centers_array, gij_array);
            },
            py::arg("domain"),
            py::arg("n_bins"),
            py::arg("maximum_pairs") = 1'000'000,
            "Compute partial g_ij(r) for this configuration in a supplied domain."
        )
        .def(
            "total_sphere_volume",
            [](const std::shared_ptr<SphereConfiguration> sphere_configuration){
                double volume = sphere_configuration->total_sphere_volume();
                py::object ureg = get_shared_ureg();
                py::object quantity = ureg.attr("Quantity")(volume, "meter ** 3");
                return quantity;
            },
            "Compute the total volume occupied by the spheres."
        )
        .def_readonly(
            "classes_index",
            &SphereConfiguration::class_index_values,
            "Integer class index for each sphere"
        )
        .def_property_readonly(
            "positions",
            [](const std::shared_ptr<SphereConfiguration> sphere_configuration) {
                py::array_t<double> output = vector3d_list_to_numpy(sphere_configuration->center_positions);
                py::object ureg = get_shared_ureg();
                py::object quantity = ureg.attr("Quantity")(output, "meter");
                return quantity;
            },
            "List of sphere center positions"
        )
        .def_property_readonly(
            "radii",
            [](const std::shared_ptr<SphereConfiguration> sphere_configuration) {
                py::array_t<double> output = double_list_to_numpy(sphere_configuration->radii_values);
                py::object ureg = get_shared_ureg();
                py::object quantity = ureg.attr("Quantity")(output, "meter");
                return quantity;
            },
            "List of sphere radii"
        )
        .def_property_readonly(
            "number_of_classes",
            [](const std::shared_ptr<SphereConfiguration> sphere_configuration) {
                if (sphere_configuration->class_index_values.empty()) return 0;

                int max_class = -1;
                for (int c : sphere_configuration->class_index_values)
                    if (c > max_class) max_class = c;
                return max_class + 1;
            },
            "Number of distinct particle radius classes"
        )
        .def("__repr__", [](const SphereConfiguration& self) {
            return "<PackingConfiguration spheres=" + std::to_string(self.radii_values.size()) + ">";
        });

    py::class_<Options, std::shared_ptr<Options>>(module, "RSAOptions", R"doc(
Stopping criteria and numerical settings for an RSA simulation.

Attributes
----------
maximum_attempts : int
    Total trial-insertion limit.
maximum_spheres : int
    Sphere-count limit; zero disables this criterion.
target_packing_fraction : float
    Target volume fraction; zero disables this criterion.
)doc")
        .def(py::init<>())
        .def_readwrite("random_seed", &Options::random_seed)
        .def_readwrite("maximum_attempts", &Options::maximum_attempts)
        .def_readwrite("maximum_spheres", &Options::maximum_spheres)
        .def_readwrite("maximum_consecutive_rejections", &Options::maximum_consecutive_rejections)
        .def_readwrite("target_packing_fraction", &Options::target_packing_fraction)
        .def_readwrite("minimum_center_separation_addition", &Options::minimum_center_separation_addition)
        .def_readwrite("containment_padding", &Options::containment_padding)
        .def_readwrite("spatial_grid_cell_size", &Options::spatial_grid_cell_size)
        .def_readwrite("enforce_radii_distribution", &Options::enforce_radii_distribution)
        .def("__repr__", [](const Options& self) {
            std::ostringstream stream;
            stream << "<RSAOptions seed=" << self.random_seed
                   << ", max_attempts=" << self.maximum_attempts
                   << ", target_packing_fraction=" << self.target_packing_fraction << ">";
            return stream.str();
        });

    py::class_<Simulator>(module, "RSASimulator", R"doc(
Random Sequential Addition simulator for non-overlapping spheres.

Parameters
----------
domain : PackingDomain
    Spatial domain and boundary conditions.
radius_sampler : RadiusSampler
    Distribution used for candidate sphere radii.
options : RSAOptions
    Stopping criteria and numerical settings.

Notes
-----
Call :meth:`run` to generate a :class:`PackingResult`.
)doc")
        .def(
            py::init<std::shared_ptr<MCDomain>, std::shared_ptr<RadiusSampler>, std::shared_ptr<Options>>(),
            py::arg("domain"),
            py::arg("radius_sampler"),
            py::arg("options")
        )
        .def(
            "reset",
            &Simulator::reset,
            "Reset the simulation to its initial state."
        )
        .def(
            "run",
            &run_and_wrap,
            "Run the simulation and return a PackingResult."
        )
        .def(
            "_cpp_attempt_single_insertion",
            &Simulator::attempt_single_insertion,
            "Attempt to insert a single sphere into the simulation."
        )
        .def_readonly(
            "_cpp_statistics",
            &Simulator::statistics,
            py::return_value_policy::reference_internal,
            "Simulation statistics"
        )
        .def_readonly(
            "sphere_configuration",
            &Simulator::sphere_configuration,
            py::return_value_policy::reference_internal,
            "Current sphere configuration"
        )
        .def("__repr__", [](const Simulator&) {
            return "<RSASimulator>";
        });


}
