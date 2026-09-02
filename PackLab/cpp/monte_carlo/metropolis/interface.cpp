#include "metropolis.h"

#include <pybind11/pybind11.h>

#include <pint/pint.h>

namespace py = pybind11;

namespace {

py::object run_and_wrap(MetropolisSimulator& simulator) {
    auto cpp_result = simulator.run();
    py::object result_class = py::module_::import("PackLab.monte_carlo.results").attr("PackingResult");
    return result_class(py::arg("binding") = py::cast(std::move(cpp_result)));
}

}  // namespace

PYBIND11_MODULE(metropolis, module) {
    module.doc() = "Equilibrium Metropolis sampling for hard-sphere configurations.";

    py::class_<MetropolisOptions, std::shared_ptr<MetropolisOptions>>(module, "MetropolisOptions", R"doc(
Numerical settings for fixed-volume hard-sphere Metropolis sampling.

The maximum displacement is the maximum absolute displacement proposed along
each coordinate in one move.  Tune it to obtain a practical acceptance rate.
)doc")
        .def(py::init<>())
        .def_readwrite("random_seed", &MetropolisOptions::random_seed)
        .def_readwrite("number_of_sweeps", &MetropolisOptions::number_of_sweeps)
        .def_property(
            "maximum_displacement",
            [](const MetropolisOptions& options) {
                return meters_quantity_with_ureg(get_shared_ureg(), options.maximum_displacement);
            },
            [](MetropolisOptions& options, const py::object& value) {
                const double displacement = to_meters_strict(value);
                if (!(displacement > 0.0) || !std::isfinite(displacement)) {
                    throw std::invalid_argument("maximum_displacement must be finite and positive.");
                }
                options.maximum_displacement = displacement;
            },
            "pint.Quantity: Maximum absolute displacement along each coordinate."
        )
        .def("__repr__", [](const MetropolisOptions& options) {
            return "<MetropolisOptions sweeps=" + std::to_string(options.number_of_sweeps) + ">";
        });

    py::class_<MetropolisStatistics>(module, "MetropolisStatistics", R"doc(
Move statistics from an equilibrium hard-sphere Metropolis simulation.

Unlike :class:`PackingStatistics`, these counters describe displacement moves,
not RSA insertion attempts.
)doc")
        .def_readonly("attempted_moves", &MetropolisStatistics::attempted_moves)
        .def_readonly("accepted_moves", &MetropolisStatistics::accepted_moves)
        .def_readonly("rejected_moves", &MetropolisStatistics::rejected_moves)
        .def_readonly("completed_sweeps", &MetropolisStatistics::completed_sweeps)
        .def_readonly("total_runtime_seconds", &MetropolisStatistics::total_runtime_seconds)
        .def_property_readonly("acceptance_rate", &MetropolisStatistics::acceptance_rate);

    py::class_<MetropolisSimulator>(module, "MetropolisSimulator", R"doc(
Equilibrate a non-overlapping hard-sphere configuration using Metropolis moves.

Parameters
----------
domain : PackingDomain
    Fixed simulation domain and boundary conditions.
initial_configuration : PackingConfiguration
    Valid, non-overlapping configuration, for example from ``RSASimulator``.
options : MetropolisOptions
    Sweep count, proposal displacement, and seed.

Notes
-----
This sampler holds particle count, radii, and class labels fixed. It is an
equilibrium workflow and is not a continuation of RSA deposition.
)doc")
        .def(
            py::init<std::shared_ptr<MCDomain>, std::shared_ptr<SphereConfiguration>, std::shared_ptr<MetropolisOptions>>(),
            py::arg("domain"),
            py::arg("initial_configuration"),
            py::arg("options")
        )
        .def("reset", &MetropolisSimulator::reset, "Restore the supplied initial configuration.")
        .def("run", &run_and_wrap, "Run all configured sweeps and return a PackingResult.")
        .def_readonly("sphere_configuration", &MetropolisSimulator::sphere_configuration,
                      py::return_value_policy::reference_internal,
                      "Current configuration after accepted displacement moves.")
        .def_readonly("statistics", &MetropolisSimulator::statistics,
                      py::return_value_policy::reference_internal,
                      "Metropolis displacement-move statistics.")
        .def("__repr__", [](const MetropolisSimulator&) { return "<MetropolisSimulator>"; });
}
