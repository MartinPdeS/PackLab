#include "estimator.h"

#include <pybind11/pybind11.h>
#include "pint/pint.h"

namespace py = pybind11;


PYBIND11_MODULE(interface_estimator, module) {

    py::class_<EstimateResult>(module, "EstimateResult")
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
