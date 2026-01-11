#include "radius_sampler.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "pint/pint.h"

namespace py = pybind11;

PYBIND11_MODULE(interface_sampler, module) {
    module.doc() = "Radius sampler classes for RSA simulations with optional radius binning.";

    py::class_<RadiusSampler, std::shared_ptr<RadiusSampler>>(module, "RadiusSampler")
        .def("set_number_of_bins", &RadiusSampler::set_number_of_bins, py::arg("bins"))
        .def_property_readonly("number_of_bins", &RadiusSampler::number_of_bins)
        .def("set_bin_spacing", &RadiusSampler::set_bin_spacing, py::arg("bin_spacing"))
        .def_property_readonly("bin_spacing", &RadiusSampler::bin_spacing)
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
        );

    py::class_<ConstantRadiusSampler, RadiusSampler, std::shared_ptr<ConstantRadiusSampler>>(module, "Constant", py::dynamic_attr())
        .def(py::init([](py::object radius_py, int bins) {
                const double radius_m = to_meters_strict(radius_py);
                return std::make_shared<ConstantRadiusSampler>(radius_m, bins);
            }),
            py::arg("radius"),
            py::arg("bins") = 0
        )
        .def_property_readonly("_ureg", [](py::object) { return get_shared_ureg(); });

    py::class_<UniformRadiusSampler, RadiusSampler, std::shared_ptr<UniformRadiusSampler>>(module, "Uniform", py::dynamic_attr())
        .def(py::init([](py::object minimum_radius_py, py::object maximum_radius_py, int bins) {
                const double minimum_radius_m = to_meters_strict(minimum_radius_py);
                const double maximum_radius_m = to_meters_strict(maximum_radius_py);
                return std::make_shared<UniformRadiusSampler>(minimum_radius_m, maximum_radius_m, bins);
            }),
            py::arg("minimum_radius"),
            py::arg("maximum_radius"),
            py::arg("bins") = 0
        )
        .def_property_readonly("_ureg", [](py::object) { return get_shared_ureg(); });

    // Important: mu and sigma are log space parameters, dimensionless.
    // If you want physical parameters, expose (median_radius, geometric_standard_deviation) instead.
    py::class_<LogNormalRadiusSampler, RadiusSampler, std::shared_ptr<LogNormalRadiusSampler>>(module, "LogNormal", py::dynamic_attr())
        .def(py::init([](double mu, double sigma, py::object maximum_radius_clip_py, int bins) {
                double maximum_radius_clip_m = 0.0;
                if (maximum_radius_clip_py.is_none()) {
                    throw py::value_error("maximum_radius_clip must be provided for this LogNormal sampler.");
                }
                maximum_radius_clip_m = to_meters_strict(maximum_radius_clip_py);

                return std::make_shared<LogNormalRadiusSampler>(mu, sigma, maximum_radius_clip_m, bins);
            }),
            py::arg("mu"),
            py::arg("sigma"),
            py::arg("maximum_radius_clip"),
            py::arg("bins") = 0
        )
        .def_property_readonly("_ureg", [](py::object) { return get_shared_ureg(); });

    py::class_<DiscreteRadiusSampler, RadiusSampler, std::shared_ptr<DiscreteRadiusSampler>>(module, "Discrete", py::dynamic_attr())
        .def(py::init([](py::object radii_py, py::object weights_py, int bins) {
                const std::vector<double> radii_m = to_vector_units(radii_py, "meter");

                py::array weights_array = py::array::ensure(weights_py);
                if (!weights_array || weights_array.ndim() != 1) {
                    throw py::value_error("weights must be a one dimensional array.");
                }
                py::array_t<double, py::array::c_style | py::array::forcecast> weights_double(weights_array);
                auto wbuf = weights_double.unchecked<1>();

                std::vector<double> weights(static_cast<std::size_t>(wbuf.shape(0)));
                for (py::ssize_t i = 0; i < wbuf.shape(0); ++i) {
                    weights[static_cast<std::size_t>(i)] = wbuf(i);
                }

                return std::make_shared<DiscreteRadiusSampler>(radii_m, weights, bins);
            }),
            py::arg("radii"),
            py::arg("weights"),
            py::arg("bins") = 0
        )
        .def_property_readonly("_ureg", [](py::object) { return get_shared_ureg(); });
}
