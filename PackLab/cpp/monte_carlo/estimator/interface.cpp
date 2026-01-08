#include "estimator.h"

#include <pybind11/pybind11.h>
#include "monte_carlo/utils/pint.h"

namespace py = pybind11;


PYBIND11_MODULE(interface_estimator, module) {

    py::class_<EstimateResult>(module, "EstimateResult")
        .def_readonly(
            "centers",
            &EstimateResult::centers
        )
        .def_readonly(
            "mean_g",
            &EstimateResult::mean_g
        )
        .def_readonly(
            "std_g",
            &EstimateResult::std_g
        )
        .def_readonly(
            "number_of_species",
            &EstimateResult::number_of_species
        )
        .def_readonly(
            "number_of_bins",
            &EstimateResult::number_of_bins
        );

    py::class_<Estimator, std::shared_ptr<Estimator>> estimator_cls(module, "Estimator");
    estimator_cls
        .def(py::init<
            const std::shared_ptr<Domain>&,
            const std::shared_ptr<RadiusSampler>&,
            const std::shared_ptr<Options>&,
            std::size_t>(),
            py::arg("domain"),
            py::arg("radius_sampler"),
            py::arg("options"),
            py::arg("number_of_bins")
        )

        .def("estimate",
            &Estimator::estimate,
            py::arg("number_of_samples"),
            py::arg("maximum_pairs") = 0
        )
    ;


}
