#include "estimator.h"

#include <pybind11/pybind11.h>
#include <sstream>
#include "pint/pint.h"

namespace py = pybind11;


PYBIND11_MODULE(estimator, module) {

    py::class_<EstimateResult>(module, "PackingEstimate", R"doc(
Ensemble estimate of partial pair-correlation functions.

Attributes
----------
centers : pint.Quantity, shape (n_bins,)
    Radial bin centers in meters.
mean_g, std_g : numpy.ndarray, shape (n_species, n_species, n_bins)
    Mean and sample standard deviation of partial pair correlations.
)doc")
        .def_property_readonly(
            "centers",
            [](const EstimateResult& self){
                py::object base = py::cast(&self, py::return_value_policy::reference);
                py::object ureg = get_shared_ureg();

                py::array_t<double> output = vector_double_to_numpy_view(self.centers, base);

                py::object quantity = ureg.attr("Quantity")(output, "meter");

                return quantity;
            }
        )
        .def_property_readonly(
            "mean_g",
            [](const EstimateResult& self){
                py::object base = py::cast(&self, py::return_value_policy::reference);
                return vector_vector_vector_double_to_numpy_plane_row_views(self.mean_g, base);
            }
        )
        .def_property_readonly(
            "std_g",
            [](const EstimateResult& self){
                py::object base = py::cast(&self, py::return_value_policy::reference);
                return vector_vector_vector_double_to_numpy_plane_row_views(self.std_g, base);
            }
        )
        .def_readonly(
            "number_of_species",
            &EstimateResult::number_of_species
        )
        .def_readonly(
            "number_of_bins",
            &EstimateResult::number_of_bins
        )
        .def("__repr__", [](const EstimateResult& self) {
            return "<PackingEstimate species=" + std::to_string(self.number_of_species)
                + ", bins=" + std::to_string(self.number_of_bins) + ">";
        });

    py::class_<Estimator, std::shared_ptr<Estimator>> estimator_cls(module, "PackingEstimator", R"doc(
Estimate pair-correlation statistics over repeated RSA simulations.

Parameters
----------
domain : PackingDomain
    Simulation domain.
radius_sampler : RadiusSampler
    Particle-radius distribution.
options : RSAOptions
    RSA simulation configuration.
number_of_bins : int
    Number of radial bins used for each estimate.
)doc");
    estimator_cls
        .def(py::init<
            const std::shared_ptr<MCDomain>&,
            const std::shared_ptr<RadiusSampler>&,
            const std::shared_ptr<Options>&,
            std::size_t>(),
            py::arg("domain"),
            py::arg("radius_sampler"),
            py::arg("options"),
            py::arg("number_of_bins"),
            "Create an ensemble packing estimator."
        )

        .def("estimate",
            &Estimator::estimate,
            py::arg("number_of_samples"),
            py::arg("maximum_pairs") = 0,
            R"doc(
Estimate partial pair correlations from independent RSA realizations.

Parameters
----------
number_of_samples : int
    Number of independent packing realizations.
maximum_pairs : int, default=0
    Pair-sampling limit per realization; zero selects the native default.

Returns
-------
PackingEstimate
    Mean and standard deviation of the partial pair-correlation functions.
)doc"
        )
        .def("__repr__", [](const Estimator&) { return "<PackingEstimator>"; })
    ;


}
