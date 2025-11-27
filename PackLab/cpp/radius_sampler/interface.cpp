#include "radius_sampler.h"

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>


PYBIND11_MODULE(interface_radius_sampler, module) {
    module.doc() = "Random sequential addition of non overlapping spheres in a 3D box";

    pybind11::class_<RadiusSampler, std::shared_ptr<RadiusSampler>>(module, "RadiusSampler");

    pybind11::class_<ConstantRadiusSampler, RadiusSampler, std::shared_ptr<ConstantRadiusSampler>>(module, "ConstantRadiusSampler")
        .def(
            pybind11::init<double>(),
            pybind11::arg("radius")
        );

    pybind11::class_<UniformRadiusSampler, RadiusSampler, std::shared_ptr<UniformRadiusSampler>>(module, "UniformRadiusSampler")
        .def(
            pybind11::init<double, double>(),
            pybind11::arg("minimum_radius"),
            pybind11::arg("maximum_radius")
        );

    pybind11::class_<LogNormalRadiusSampler, RadiusSampler, std::shared_ptr<LogNormalRadiusSampler>>(module, "LogNormalRadiusSampler")
        .def(
            pybind11::init<double, double, double>(),
            pybind11::arg("mu"),
            pybind11::arg("sigma"),
            pybind11::arg("maximum_radius_clip")
        );

    pybind11::class_<DiscreteRadiusSampler, RadiusSampler, std::shared_ptr<DiscreteRadiusSampler>>(module, "DiscreteRadiusSampler")
        .def(
            pybind11::init<std::vector<double>, std::vector<double>>(),
            pybind11::arg("radii"),
            pybind11::arg("weights")
        )
        ;

    pybind11::class_<PythonCallableRadiusSampler, RadiusSampler, std::shared_ptr<PythonCallableRadiusSampler>>(module, "PythonCallableRadiusSampler")
        .def(
            pybind11::init<std::function<double()>, double>(),
            pybind11::arg("python_callable"),
            pybind11::arg("maximum_possible_radius")
        );

}
