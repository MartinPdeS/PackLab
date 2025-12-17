#include "radius_sampler.h"

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

PYBIND11_MODULE(interface_radius_sampler, module) {
    module.doc() = "Radius sampler classes for RSA simulations with optional radius binning.";

    // Base class -----------------------------------------------------
    pybind11::class_<RadiusSampler, std::shared_ptr<RadiusSampler>>(module, "RadiusSampler")
        .def(
            "set_number_of_bins",
            &RadiusSampler::set_number_of_bins,
            pybind11::arg("bins"),
            "Enable optional radius binning. bins=0 disables binning."
        )
        .def_property_readonly(
            "number_of_bins",
            &RadiusSampler::number_of_bins,
            "Return the number of bins used (0 means no binning)."
        );

    // Constant sampler ------------------------------------------------
    pybind11::class_<ConstantRadiusSampler, RadiusSampler, std::shared_ptr<ConstantRadiusSampler>>(
        module,
        "ConstantRadiusSampler"
    )
        .def(
            pybind11::init<double, int>(),
            pybind11::arg("radius"),
            pybind11::arg("bins") = 0,
            "Constant radius sampler with optional binning."
        );

    // Uniform sampler -------------------------------------------------
    pybind11::class_<UniformRadiusSampler, RadiusSampler, std::shared_ptr<UniformRadiusSampler>>(
        module,
        "UniformRadiusSampler"
    )
        .def(
            pybind11::init<double, double, int>(),
            pybind11::arg("minimum_radius"),
            pybind11::arg("maximum_radius"),
            pybind11::arg("bins") = 0,
            "Uniform radius sampler with optional binning."
        );

    // Lognormal sampler -----------------------------------------------
    pybind11::class_<LogNormalRadiusSampler, RadiusSampler, std::shared_ptr<LogNormalRadiusSampler>>(
        module,
        "LogNormalRadiusSampler"
    )
        .def(
            pybind11::init<double, double, double, int>(),
            pybind11::arg("mu"),
            pybind11::arg("sigma"),
            pybind11::arg("maximum_radius_clip"),
            pybind11::arg("bins") = 0,
            "Log-normal radius sampler with clipping and optional binning."
        );

    // Discrete sampler ------------------------------------------------
    pybind11::class_<DiscreteRadiusSampler, RadiusSampler, std::shared_ptr<DiscreteRadiusSampler>>(
        module,
        "DiscreteRadiusSampler"
    )
        .def(
            pybind11::init<std::vector<double>, std::vector<double>>(),
            pybind11::arg("radii"),
            pybind11::arg("weights"),
            "Discrete radius sampler with user-provided radii and weights plus optional binning."
        );
}
