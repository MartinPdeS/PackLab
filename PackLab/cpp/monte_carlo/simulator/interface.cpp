#include "simulator.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

py::object run_and_wrap(Simulator& self) {
    auto cpp_result = self.run();  // whatever your C++ returns

    py::object ResultClass = py::module_::import("PackLab.monte_carlo.results").attr("Result");

    return ResultClass(py::arg("binding") = py::cast(std::move(cpp_result)));
}

PYBIND11_MODULE(interface_simulator, module) {
    module.doc() = "Random sequential addition of non overlapping spheres in a 3D box";

    py::class_<SphereConfiguration, std::shared_ptr<SphereConfiguration>>(module, "SphereConfiguration")
        .def_property_readonly(
            "count",
            [](const std::shared_ptr<SphereConfiguration> sphere_configuration) { return sphere_configuration->radii().size(); },
            "Number of spheres in the configuration."
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
    ;

    py::class_<Options, std::shared_ptr<Options>>(module, "Options")
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
    ;

    py::class_<Simulator>(module, "Simulator")
        .def(
            py::init<std::shared_ptr<Domain>, std::shared_ptr<RadiusSampler>, std::shared_ptr<Options>>(),
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
            &run_and_wrap
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
        ;


}

