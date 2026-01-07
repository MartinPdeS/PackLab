#include "radius_sampler.h"

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

#include "utils/pint.h"

namespace py = pybind11;

PYBIND11_MODULE(interface_radius_sampler, module) {
    module.doc() = "Radius sampler classes for RSA simulations with optional radius binning.";

    // Base class -----------------------------------------------------
    py::class_<RadiusSampler, std::shared_ptr<RadiusSampler>>(module, "RadiusSampler")
        .def(
            "set_number_of_bins",
            &RadiusSampler::set_number_of_bins,
            py::arg("bins"),
            "Enable optional radius binning. bins=0 disables binning."
        )
        .def_property_readonly(
            "number_of_bins",
            &RadiusSampler::number_of_bins,
            "Return the number of bins used (0 means no binning)."
        );

    // Constant sampler ------------------------------------------------
    py::class_<ConstantRadiusSampler, RadiusSampler, std::shared_ptr<ConstantRadiusSampler>>(module, "Constant", py::dynamic_attr())
        .def(
            "__init__",
            [](
                py::object& self,
                py::object radius_py,
                int bins
            ) {
                py::object ureg = registry_from_quantity(radius_py);
                const double radius = to_meters_strict(radius_py);

                new (self.cast<ConstantRadiusSampler*>()) ConstantRadiusSampler(radius, bins);
                self.attr("_ureg") = ureg;
            }
        )
        .def(
            py::init<double, int>(),
            py::arg("radius"),
            py::arg("bins") = 0,
            "Constant radius sampler with optional binning."
        );

    // Uniform sampler -------------------------------------------------
    py::class_<UniformRadiusSampler, RadiusSampler, std::shared_ptr<UniformRadiusSampler>>(module, "Uniform", py::dynamic_attr())
        .def(
            "__init__",
            [](
                py::object& self,
                py::object minimum_radius_py,
                py::object maximum_radius_py,
                int bins
            ) {
                py::object ureg = registry_from_quantity(minimum_radius_py);

                const double minimum_radius = to_meters_strict(minimum_radius_py);
                const double maximum_radius = to_meters_strict(maximum_radius_py);

                new (self.cast<UniformRadiusSampler*>()) UniformRadiusSampler(minimum_radius, maximum_radius, bins);
                self.attr("_ureg") = ureg;
            },
            py::arg("minimum_radius"),
            py::arg("maximum_radius"),
            py::arg("bins") = 0,
            "Uniform radius sampler with optional binning."
        );

    // Lognormal sampler -----------------------------------------------
    py::class_<LogNormalRadiusSampler, RadiusSampler, std::shared_ptr<LogNormalRadiusSampler>>(module, "LogNormal", py::dynamic_attr())
        .def(
            "__init__",
            [](
                py::object self,
                py::object mu_py,
                py::object sigma_py,
                py::object maximum_radius_clip_py,
                int bins
            ) {
                py::object ureg = registry_from_quantity(mu_py);

                const double mu = to_meters_strict(mu_py);
                const double sigma = to_meters_strict(sigma_py);

                double maximum_radius_clip = 0.0;
                if (maximum_radius_clip_py.is_none()) {
                    // Example default: 10 * mu (choose what makes sense for your model)
                    maximum_radius_clip = 10.0 * mu;
                } else {
                    maximum_radius_clip = to_meters_strict(maximum_radius_clip_py);
                }

                new (self.cast<LogNormalRadiusSampler*>()) LogNormalRadiusSampler(mu, sigma, maximum_radius_clip, bins);
                self.attr("_ureg") = ureg;
            },
            py::arg("mu"),
            py::arg("sigma"),
            py::arg("maximum_radius_clip") = py::none(),
            py::arg("bins") = 0,
            "Log-normal radius sampler with clipping and optional binning."
        );


    // Discrete sampler ------------------------------------------------
    py::class_<DiscreteRadiusSampler, RadiusSampler, std::shared_ptr<DiscreteRadiusSampler>>(module, "Discrete", py::dynamic_attr())
        .def(
            "__init__",
            [](
                py::object& self,
                py::object radii_py,
                std::vector<double> weights
            ) {
                py::object ureg = registry_from_quantity(radii_py);
                const std::vector<double> radii = to_vector_units(radii_py, "meter");
                new (self.cast<DiscreteRadiusSampler*>()) DiscreteRadiusSampler(radii, weights);
                self.attr("_ureg") = ureg;
            },
            py::arg("radii"),
            py::arg("weights"),
            "Discrete radius sampler with user-provided radii and weights plus optional binning."
        )
        ;
}
