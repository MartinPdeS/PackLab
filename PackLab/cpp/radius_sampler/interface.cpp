#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <sstream>

#include "pint/pint.h"
#include "uniform.h"
#include "normal.h"
#include "log_normal.h"
#include "discrete.h"
#include "constant.h"

namespace py = pybind11;

PYBIND11_MODULE(samplers, module) {
    py::object ureg = get_shared_ureg();

    module.doc() = "Radius sampler classes for RSA simulations with optional radius binning.";



    // Base class -----------------------------------------------------
    py::class_<RadiusSampler, std::shared_ptr<RadiusSampler>>(module, "RadiusSampler", R"doc(
Base class for particle-radius distributions.

Notes
-----
Concrete samplers generate radii in SI units internally. Their Python
constructors accept Pint quantities and ``to_bins`` returns radii in meters.
)doc")
        .def(
            "set_number_of_bins",
            &RadiusSampler::set_number_of_bins,
            py::arg("bins"),
            R"doc(
Set the number of discrete radius bins.

Parameters
----------
bins : int
    Number of bins. A value of zero disables binning.
)doc"
        )
        .def_property_readonly(
            "number_of_bins",
            &RadiusSampler::number_of_bins,
            "int: Number of bins; zero means that binning is disabled."
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
            },
            R"doc(
Convert the distribution to representative discrete radius classes.

Returns
-------
radii : pint.Quantity, shape (n_bins,)
    Representative radii in meters.
weights : numpy.ndarray, shape (n_bins,)
Normalized probability weights.
)doc"
        )
        .def(
            "plot_histogram",
            [ureg](const RadiusSampler& self, py::object ax, py::object unit, py::kwargs kwargs) {
                auto [radii_m, weights] = self.to_bins();

                py::array_t<double> radii_array(static_cast<py::ssize_t>(radii_m.size()));
                py::array_t<double> weights_array(static_cast<py::ssize_t>(weights.size()));
                {
                    auto radii_buffer = radii_array.mutable_unchecked<1>();
                    auto weights_buffer = weights_array.mutable_unchecked<1>();
                    for (py::ssize_t i = 0; i < radii_buffer.shape(0); ++i) {
                        radii_buffer(i) = radii_m[static_cast<std::size_t>(i)];
                        weights_buffer(i) = weights[static_cast<std::size_t>(i)];
                    }
                }

                py::object radii = (radii_array * ureg.attr("meter")).attr("to")(unit).attr("magnitude");
                if (ax.is_none()) {
                    ax = py::module_::import("matplotlib.pyplot").attr("gca")();
                }

                py::tuple arguments(2);
                arguments[0] = radii;
                arguments[1] = weights_array;
                ax.attr("bar")(*arguments, **kwargs);
                ax.attr("set_xlabel")(py::str("particle radius (") + py::str(unit) + py::str(")"));
                ax.attr("set_ylabel")("number fraction");
                return ax;
            },
            py::arg("ax") = py::none(),
            py::arg("unit") = "nanometer",
            R"doc(
Plot the binned radius distribution as a probability-mass histogram.

The bars show the representative radius classes and normalized number
fractions returned by :meth:`to_bins`. Continuous samplers therefore require
``bins`` to be positive.

Parameters
----------
ax : matplotlib.axes.Axes, optional
    Axes receiving the bars. The current axes are used when omitted.
unit : str or pint.Unit, default="nanometer"
    Unit used on the radius axis.
**kwargs
    Additional keyword arguments forwarded to :meth:`matplotlib.axes.Axes.bar`.

Returns
-------
matplotlib.axes.Axes
    The axes containing the histogram.
)doc"
        )
        ;

    // Constant sampler ------------------------------------------------
    py::class_<Constant, RadiusSampler, std::shared_ptr<Constant>>(module, "ConstantRadiusSampler", R"doc(
Radius sampler that always returns one radius.

Parameters
----------
radius : pint.Quantity
    Particle radius.
bins : int, default=0
    Optional number of radius bins.
)doc")
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
            "Create a constant-radius sampler."
        )
        .def("__repr__", [](const Constant& self) {
            return "<ConstantRadiusSampler bins=" + std::to_string(self.number_of_bins()) + ">";
        })
        ;


    // Uniform sampler -------------------------------------------------
    py::class_<Uniform, RadiusSampler, std::shared_ptr<Uniform>>(module, "UniformRadiusSampler", R"doc(
Uniform distribution of particle radii.

Parameters
----------
minimum_radius, maximum_radius : pint.Quantity
    Inclusive bounds of the radius distribution.
bins : int, default=0
    Optional number of radius bins.
)doc")
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
            "Create a uniform radius sampler."
        )
        .def("__repr__", [](const Uniform& self) {
            return "<UniformRadiusSampler bins=" + std::to_string(self.number_of_bins()) + ">";
        })
        ;

    // Lognormal sampler -----------------------------------------------
    py::class_<LogNormal, RadiusSampler, std::shared_ptr<LogNormal>>(module, "LogNormalRadiusSampler", R"doc(
Log-normal distribution of strictly positive particle radii.

Parameters
----------
median_radius : pint.Quantity
    Median radius.
geometric_standard_deviation : float
    Multiplicative standard deviation; must exceed one.
maximum_radius_clip : pint.Quantity
    Upper radius cutoff.
bins : int, default=0
    Optional number of radius bins.
)doc")
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
        .def("__repr__", [](const LogNormal& self) {
            return "<LogNormalRadiusSampler bins=" + std::to_string(self.number_of_bins()) + ">";
        })
        ;



    // Discrete sampler ------------------------------------------------
    py::class_<Discrete, RadiusSampler, std::shared_ptr<Discrete>>(module, "DiscreteRadiusSampler", R"doc(
Finite set of particle radii with explicit probability weights.

Parameters
----------
radii : pint.Quantity, shape (n,)
    Discrete radius values.
weights : array-like, shape (n,)
    Non-negative relative probabilities. They are normalized internally.
)doc")
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
            "Create a discrete radius sampler."
        )
        .def("__repr__", [](const Discrete& self) {
            return "<DiscreteRadiusSampler bins=" + std::to_string(self.number_of_bins()) + ">";
        })
        ;

    // Normal (Gaussian) sampler ----------------------------------------------
    py::class_<Normal, RadiusSampler, std::shared_ptr<Normal>>(module, "NormalRadiusSampler", R"doc(
Clipped normal distribution of particle radii.

Parameters
----------
mean, standard_deviation : pint.Quantity
    Parameters of the normal distribution.
maximum_clip : pint.Quantity, optional
    Upper radius cutoff. Defaults to ``mean + 5 * standard_deviation``.
bins : int, default=0
    Optional number of radius bins.
)doc")
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
                    }
                ),
                py::arg("mean"),
                py::arg("standard_deviation"),
                py::arg("maximum_clip") = py::none(),
                py::arg("bins") = 0,
                "Create a clipped normal radius sampler."
            )
            .def("__repr__", [](const Normal& self) {
                return "<NormalRadiusSampler bins=" + std::to_string(self.number_of_bins()) + ">";
            })
        ;

}
