#include "rsa.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

static pybind11::array_t<double> vector3d_list_to_numpy(const std::vector<Vector3d>& values) {
    const pybind11::ssize_t count = static_cast<pybind11::ssize_t>(values.size());
    auto array = pybind11::array_t<double>({count, pybind11::ssize_t(3)});
    auto buffer = array.mutable_unchecked<2>();

    for (pybind11::ssize_t index = 0; index < count; ++index) {
        buffer(index, 0) = values[static_cast<std::size_t>(index)].x;
        buffer(index, 1) = values[static_cast<std::size_t>(index)].y;
        buffer(index, 2) = values[static_cast<std::size_t>(index)].z;
    }

    return array;
}

static pybind11::array_t<double> double_list_to_numpy(const std::vector<double>& values) {
    const pybind11::ssize_t count = static_cast<pybind11::ssize_t>(values.size());
    auto array = pybind11::array_t<double>(count);
    std::memcpy(array.mutable_data(), values.data(), static_cast<std::size_t>(count) * sizeof(double));
    return array;
}

PYBIND11_MODULE(interface_rsa, module) {
    module.doc() = "Random sequential addition of non overlapping spheres in a 3D box";

    pybind11::class_<SphereConfiguration>(module, "SphereConfiguration")
        .def_property_readonly(
            "count",
            [](const SphereConfiguration& configuration) { return configuration.radii().size(); },
            "Number of spheres in the configuration."
        )
        .def(
            "total_sphere_volume",
            &SphereConfiguration::total_sphere_volume,
            "Compute the total volume occupied by the spheres."
        )
        .def_readonly(
            "classes_index",
            &SphereConfiguration::class_index_values,
            "Integer class index for each sphere"
        )
        .def_property_readonly(
            "positions",
            [](const SphereConfiguration& sphere_configuration) {
                return vector3d_list_to_numpy(sphere_configuration.center_positions);
            },
            "List of sphere center positions"
        )
        .def_property_readonly(
            "radii",
            [](const SphereConfiguration& sphere_configuration) {
                return double_list_to_numpy(sphere_configuration.radii_values);
            },
            "List of sphere radii"
        )
        .def_property_readonly(
            "number_of_classes",
            [](const SphereConfiguration& config) {
                if (config.class_index_values.empty()) return 0;

                int max_class = -1;
                for (int c : config.class_index_values)
                    if (c > max_class) max_class = c;
                return max_class + 1;
            },
            "Number of distinct particle radius classes"
        )
    ;

    pybind11::class_<Options>(module, "Options")
        .def(pybind11::init<>())
        .def_readwrite("random_seed", &Options::random_seed)
        .def_readwrite("maximum_attempts", &Options::maximum_attempts)
        .def_readwrite("maximum_spheres", &Options::maximum_spheres)
        .def_readwrite("maximum_consecutive_rejections", &Options::maximum_consecutive_rejections)
        .def_readwrite("target_packing_fraction", &Options::target_packing_fraction)
        .def_readwrite("minimum_center_separation_addition", &Options::minimum_center_separation_addition)
        .def_readwrite("containment_padding", &Options::containment_padding)
        .def_readwrite("spatial_grid_cell_size", &Options::spatial_grid_cell_size)
        .def_readwrite("enforce_radii_distribution", &Options::enforce_radii_distribution)
    ;

    pybind11::class_<Simulator>(module, "Simulator")
        .def(
            pybind11::init<Domain, std::shared_ptr<RadiusSampler>, Options>(),
            pybind11::arg("domain"),
            pybind11::arg("radius_sampler"),
            pybind11::arg("options")
        )
        .def(
            "_cpp_reset",
            &Simulator::reset,
            "Reset the simulation to its initial state."
        )
        .def(
            "_cpp_run",
            &Simulator::run
        )
        .def(
            "_cpp_attempt_single_insertion",
            &Simulator::attempt_single_insertion,
            "Attempt to insert a single sphere into the simulation."
        )
        .def_readonly(
            "_cpp_statistics",
            &Simulator::statistics,
            pybind11::return_value_policy::reference_internal,
            "Simulation statistics"
        )
        .def_readonly(
            "domain",
            &Simulator::domain,
            pybind11::return_value_policy::reference_internal,
            "Simulation domain"
        )
        .def_readonly(
            "sphere_configuration",
            &Simulator::sphere_configuration,
            pybind11::return_value_policy::reference_internal,
            "Current sphere configuration"
        )
        ;


}

