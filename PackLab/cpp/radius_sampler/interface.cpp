#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

#include "pint/pint.h"
#include "uniform.h"
#include "normal.h"
#include "log_normal.h"
#include "discrete.h"
#include "constant.h"

namespace py = pybind11;

PYBIND11_MODULE(interface_radius_sampler, module) {
    py::object ureg = get_shared_ureg();

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
            [ureg](const RadiusSampler& self) {


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
    py::class_<Constant, RadiusSampler, std::shared_ptr<Constant>>(module, "Constant")
        .def(
            py::init(
                [ureg](
                    py::object radius_py,
                    int bins
                ) {
                    const double radius = to_meters_strict(radius_py);
                    auto sampler = std::make_shared<Constant>(radius, bins);
                    return sampler;
                }
            ),
            py::arg("radius"),
            py::arg("bins") = 0,
            "Constant radius sampler with optional binning."
        )
        ;


    // Uniform sampler -------------------------------------------------
    py::class_<Uniform, RadiusSampler, std::shared_ptr<Uniform>>(module, "Uniform")
        .def(
            py::init(
                [ureg](
                    py::object minimum_radius_py,
                    py::object maximum_radius_py,
                    int bins
                ) {

                    const double minimum_radius = to_meters_strict(minimum_radius_py);
                    const double maximum_radius = to_meters_strict(maximum_radius_py);

                    return std::make_shared<Uniform>(minimum_radius, maximum_radius, bins);
                }
            ),
            py::arg("minimum_radius"),
            py::arg("maximum_radius"),
            py::arg("bins") = 0,
            "Uniform radius sampler with optional binning."
        )
        ;

    // Lognormal sampler -----------------------------------------------
    py::class_<LogNormal, RadiusSampler, std::shared_ptr<LogNormal>>(module, "LogNormal")
        .def(
            py::init(
                [ureg](
                    py::object median_radius_py,
                    double geometric_standard_deviation,
                    py::object maximum_radius_clip_py,
                    int bins
                ) {

                    // const double mediam_radius_m = to_meters_strict(median_radius_py);
                    const double median_radius_m = median_radius_py.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                    if (!(geometric_standard_deviation > 1.0)) {
                        throw py::value_error("geometric_standard_deviation must be > 1.");
                    }

                    const double mu = std::log(median_radius_m);
                    const double sigma = std::log(geometric_standard_deviation);

                    const double maximum_radius_clip_m = to_meters_strict(maximum_radius_clip_py);

                    return std::make_shared<LogNormal>(mu, sigma, maximum_radius_clip_m, bins);
                }
            ),
            py::arg("median_radius"),
            py::arg("geometric_standard_deviation"),
            py::arg("maximum_radius_clip"),
            py::arg("bins") = 0
        )
        ;



    // Discrete sampler ------------------------------------------------
    py::class_<Discrete, RadiusSampler, std::shared_ptr<Discrete>>(module, "Discrete")
        .def(
            py::init(
                [ureg](
                    py::object radii_py,
                    std::vector<double> weights
                ) {
                    const std::vector<double> radii = to_vector_units(radii_py, "meter");

                    if (radii.size() != weights.size()) {
                        throw py::value_error("radii and weights must have the same length.");
                    }

                    return std::make_shared<Discrete>(radii, weights);
                }
            ),
            py::arg("radii"),
            py::arg("weights"),
            "Discrete radius sampler with user-provided radii and weights plus optional binning."
        )
        ;

    // Normal (Gaussian) sampler ----------------------------------------------
    py::class_<Normal, RadiusSampler, std::shared_ptr<Normal>>(module, "Normal")
            .def(
                py::init(
                    [ureg](
                        py::object mean_py,
                        py::object standard_deviation_py,
                        py::object maximum_clip_py,
                        int bins
                    ) {
                        const double mean_radius_m = mean_py.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();
                        const double sigma_radius_m = standard_deviation_py.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                        if (!std::isfinite(sigma_radius_m) || sigma_radius_m < 0.0) {
                            throw py::value_error("standard_deviation must be >= 0 and finite.");
                        }

                        if (maximum_clip_py.is(py::none()))
                            maximum_clip_py = mean_py + py::float_(5.0) * standard_deviation_py;

                        const double maximum_clip_m = maximum_clip_py.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();
                        return std::make_shared<Normal>(mean_radius_m, sigma_radius_m, maximum_clip_m, bins);
                        return std::make_shared<Normal>(1, 1, 1, 1);
                    }
                ),
                py::arg("mean"),
                py::arg("standard_deviation"),
                py::arg("maximum_clip") = py::none(),
                py::arg("bins") = 0,
                "Gaussian (normal) radius sampler. Samples N(mean, standard_deviation) in meters, then clips to [0, maximum_clip], with optional binning."
            )
        ;

}
