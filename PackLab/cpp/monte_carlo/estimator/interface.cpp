#include "estimator.h"

#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <sstream>
#include "pint/pint.h"

namespace py = pybind11;


PYBIND11_MODULE(estimator, module) {

    py::class_<EstimatorStatistics>(module, "EstimatorStatistics", R"doc(
Aggregate diagnostics from the most recent ``PackingEstimator.estimate`` call.

Attributes
----------
requested_samples, completed_samples : int
    Requested and successfully completed RSA realisations.
attempted_insertions, accepted_insertions, rejected_insertions : int
    Totals over all completed realisations.
acceptance_rate : float
    Accepted insertions divided by attempted insertions.
mean_sphere_count, mean_packing_fraction : float
    Per-realisation means.
total_runtime_seconds, mean_runtime_seconds : float
    Wall-clock simulation times.
)doc")
        .def_readonly("requested_samples", &EstimatorStatistics::requested_samples)
        .def_readonly("completed_samples", &EstimatorStatistics::completed_samples)
        .def_readonly("attempted_insertions", &EstimatorStatistics::attempted_insertions)
        .def_readonly("accepted_insertions", &EstimatorStatistics::accepted_insertions)
        .def_readonly("rejected_insertions", &EstimatorStatistics::rejected_insertions)
        .def_readonly("total_spheres", &EstimatorStatistics::total_spheres)
        .def_readonly("total_runtime_seconds", &EstimatorStatistics::total_runtime_seconds)
        .def_readonly("mean_packing_fraction", &EstimatorStatistics::mean_packing_fraction)
        .def_property_readonly("acceptance_rate", &EstimatorStatistics::acceptance_rate)
        .def_property_readonly("mean_sphere_count", &EstimatorStatistics::mean_sphere_count)
        .def_property_readonly("mean_runtime_seconds", &EstimatorStatistics::mean_runtime_seconds)
        .def("print", [](const EstimatorStatistics& self) {
            py::scoped_ostream_redirect redirect(std::cout, py::module_::import("sys").attr("stdout"));
            self.print();
        }, "Print a tabular summary of the aggregate diagnostics.")
        .def("__repr__", [](const EstimatorStatistics& self) {
            return "<EstimatorStatistics completed_samples=" + std::to_string(self.completed_samples)
                + "/" + std::to_string(self.requested_samples) + ">";
        });

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

Attributes
----------
statistics : EstimatorStatistics
    Aggregate diagnostics from the most recent call to ``estimate``.
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
            [](Estimator& self,
               const std::size_t number_of_samples,
               const std::size_t maximum_pairs,
               const bool progress,
               const std::size_t progress_interval) {
                // Route native progress output through ``sys.stdout``.  In
                // particular, this makes pytest's Windows capture mechanism
                // observe output from C++ ``std::cout``.
                py::scoped_ostream_redirect redirect(std::cout, py::module_::import("sys").attr("stdout"));
                return self.estimate(number_of_samples, maximum_pairs, progress, progress_interval);
            },
            py::arg("number_of_samples"),
            py::arg("maximum_pairs") = 0,
            py::arg("progress") = false,
            py::arg("progress_interval") = 1,
            R"doc(
Estimate partial pair correlations from independent RSA realizations.

Parameters
----------
number_of_samples : int
    Number of independent packing realizations.
maximum_pairs : int, default=0
    Pair-sampling limit per realization; zero selects the native default.
progress : bool, default=False
    Print completed-sample progress and per-sample insertion diagnostics.
progress_interval : int, default=1
    Print every this many completed samples when ``progress=True``.

Returns
-------
PackingEstimate
    Mean and standard deviation of the partial pair-correlation functions.
)doc"
        )
        .def_property_readonly(
            "statistics",
            [](Estimator& self) -> EstimatorStatistics& { return self.statistics; },
            py::return_value_policy::reference_internal,
            "Aggregate diagnostics from the most recent estimate."
        )
        .def("print_statistics", [](const Estimator& self) {
            py::scoped_ostream_redirect redirect(std::cout, py::module_::import("sys").attr("stdout"));
            self.statistics.print();
        },
            "Print aggregate diagnostics from the most recent estimate.")
        .def("__repr__", [](const Estimator&) { return "<PackingEstimator>"; })
    ;


}
