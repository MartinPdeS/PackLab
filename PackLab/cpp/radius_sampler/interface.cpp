#include "radius_sampler.h"

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

#include "pint/pint.h"

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
        )
        .def(
            "to_bins",
            [](const RadiusSampler& self) {
                py::object ureg = get_shared_ureg();

                auto [radii_m, weights] = self.to_bins();

                py::array_t<double> radii_array(static_cast<py::ssize_t>(radii_m.size()));
                {
                    auto buf = radii_array.mutable_unchecked<1>();
                    for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                        buf(i) = radii_m[static_cast<std::size_t>(i)];
                    }
                }

                py::array_t<double> weights_array(static_cast<py::ssize_t>(weights.size()));
                {
                    auto buf = weights_array.mutable_unchecked<1>();
                    for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                        buf(i) = weights[static_cast<std::size_t>(i)];
                    }
                }

                py::object radii_quantity = radii_array * ureg.attr("meter");
                return py::make_tuple(radii_quantity, weights_array);
            }
        )
        ;

    // Constant sampler ------------------------------------------------
    py::class_<ConstantRadiusSampler, RadiusSampler, std::shared_ptr<ConstantRadiusSampler>>(module, "Constant", py::dynamic_attr())
        .def(
            "__init__",
            [](
                py::object& self,
                py::object radius_py,
                int bins
            ) {
                py::object ureg = get_shared_ureg();
                const double radius = to_meters_strict(radius_py);

                new (self.cast<ConstantRadiusSampler*>()) ConstantRadiusSampler(radius, bins);
                self.attr("_ureg") = ureg;
            },
            py::arg("radius"),
            py::arg("bins") = 0
        )
        ;


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
                py::object ureg = get_shared_ureg();

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
            py::init([](py::object median_radius_py,
                        double geometric_standard_deviation,
                        py::object maximum_radius_clip_py,
                        int bins) {

                const double median_radius_m = to_meters_strict(median_radius_py);

                if (!(geometric_standard_deviation > 1.0)) {
                    throw py::value_error("geometric_standard_deviation must be > 1.");
                }

                const double mu = std::log(median_radius_m);
                const double sigma = std::log(geometric_standard_deviation);

                const double maximum_radius_clip_m = to_meters_strict(maximum_radius_clip_py);

                return std::make_shared<LogNormalRadiusSampler>(mu, sigma, maximum_radius_clip_m, bins);
            }),
            py::arg("median_radius"),
            py::arg("geometric_standard_deviation"),
            py::arg("maximum_radius_clip"),
            py::arg("bins") = 0
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
                py::object ureg = get_shared_ureg();
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
