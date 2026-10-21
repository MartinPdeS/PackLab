#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <sstream>

#include "analytical/domain_analytical/domain_analytical.h"
#include "pint/pint.h"

namespace py = pybind11;

PYBIND11_MODULE(domain, module) {
    module.doc() = R"doc(
            Polydisperse cubic domain (analytical) with unit handling performed in the pybind11 wrapper.

            This module exposes a C++ implementation of a cubic simulation domain containing a
            polydisperse population of spherical particles. The numerical core operates in
            base SI units (meters), while the Python interface accepts and returns Pint quantities.

            Key conventions
            ---------------
            - All lengths are stored internally in meters.
            - Radii arrays are returned as Pint quantities in meters.
            - Volumes are returned as Pint quantities in meter**3.
            - Number densities are returned as Pint quantities in 1/meter**3.
            - Fractions (number fractions and volume fractions) are dimensionless NumPy arrays.

            The domain is specified by:
            - size: cubic side length
            - radii: per bin particle radii
            - volume_fraction: total occupied volume fraction
            - number_fractions: per bin number fractions (normalized to sum to 1)

            The class provides deterministic conversion from fractional counts to integer
            particle counts using the configured rounding mode.
        )doc"
        ;

        py::enum_<PYDomain::RoundingMode>(
            module,
            "RoundingMode",
            R"doc(
                Rounding mode used when converting expected (non integer) particle counts to integers.

                ``floor`` always rounds down. ``round`` rounds to the nearest
                integer. The selected mode affects particle counts and densities
                derived from the supplied number fractions.
            )doc"
        )
        .value(
            "floor",
            PYDomain::RoundingMode::
            Floor
        )
        .value(
            "round",
            PYDomain::RoundingMode::
            Round
        );

    py::class_<PYDomain, std::shared_ptr<PYDomain>>(
        module,
        "PercusYevickDomain",
        py::dynamic_attr(),
        R"doc(
            Cubic simulation domain containing a polydisperse population of spherical particles.

            The numerical core stores all lengths in meters and all fractions as dimensionless
            floating point values. The wrapper converts Pint quantities to meters on input and
            reconstructs Pint quantities on output.

            Parameters
            ----------
            size : pint.Quantity
                Side length of the cubic domain. Converted to meters internally.
            radii : pint.Quantity[ndarray]
                One dimensional array of particle radii defining the bins. Converted to meters internally.
            volume_fraction : float
                Total occupied volume fraction (dimensionless). Typical range is (0, 1].
            number_fractions : array_like of float
                Number fractions per bin. Must be non negative and sum to a positive value. The
                constructor normalizes them to sum to 1.
            rounding_mode : RoundingMode, default=floor
                Rounding mode used when converting expected counts to integers.

            Notes
            -----
            The model uses number fractions and a total volume fraction as its
            sole specification mechanism. Alternative specifications should use
            separate constructors. Returned arrays are NumPy arrays; quantities
            use SI units (meter, meter**3, 1/meter**3).
            )doc"
        )
        .def(
            py::init([](
                py::object size_py,
                py::object radii_py,
                py::object volume_fraction_py,
                py::object number_fractions_py,
                PYDomain::RoundingMode rounding_mode
            ) {
                const double size_m = quantity_scalar_to_meters(size_py);
                const std::vector<double> radii_m = quantity_1d_to_meters_vector(radii_py);

                const double volume_fraction = py::float_(volume_fraction_py);
                const std::vector<double> number_fraction = array_like_1d_to_double_vector(number_fractions_py);

                return PYDomain(size_m, radii_m, volume_fraction, number_fraction, rounding_mode);
            }),
            py::arg("size"),
            py::arg("radii"),
            py::arg("volume_fraction"),
            py::arg("number_fractions"),
            py::arg("rounding_mode") = PYDomain::RoundingMode::Floor,
            R"doc(
                Create a polydisperse cubic domain.

                All unit conversion is performed by the wrapper:
                - `size` and `radii` must be Pint quantities convertible to meters.
                - `volume_fraction` and `number_fractions` are dimensionless.

                Returns
                -------
                PYDomain
                    A configured domain instance.
            )doc"
        )
        .def_property_readonly(
            "size",
            [](const PYDomain& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.size) * ureg.attr("meter");
            },
            R"doc(
                PYDomain side length.

                Returns
                -------
                pint.Quantity
                    The cubic side length as a Pint quantity in meters.
            )doc"
        )
        .def_property_readonly(
            "radii",
            [](const PYDomain& self) {
                py::array_t<double> arr(static_cast<py::ssize_t>(self.radii.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = self.radii[static_cast<std::size_t>(i)];
                }

                py::object ureg = get_shared_ureg();
                return arr * ureg.attr("meter");
            },
            R"doc(
                Particle radii bins.

                Returns
                -------
                pint.Quantity[ndarray]
                    One dimensional Pint quantity array of radii in meters.
            )doc"
        )
        .def_readonly(
            "volume_fraction",
            &PYDomain::volume_fraction,
            R"doc(
                Total occupied volume fraction.

                Returns
                -------
                float
                    Dimensionless occupied volume fraction used to infer total particle volume.
            )doc"
        )
        .def_property_readonly(
            "number_fractions",
            [](const PYDomain& self) {
                py::array_t<double> arr(static_cast<py::ssize_t>(self.number_fractions.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = self.number_fractions[static_cast<std::size_t>(i)];
                }
                return arr;
            },
            R"doc(
                Per bin number fractions.

                Returns
                -------
                numpy.ndarray
                    One dimensional array of number fractions (dimensionless). Always normalized to sum to 1.
            )doc"
        )
        .def_property_readonly(
            "volume",
            [](const PYDomain& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_volume()) * ureg.attr("meter**3");
            },
            R"doc(
                PYDomain volume.

                Returns
                -------
                pint.Quantity
                    Total cubic volume in meter**3.
            )doc"
        )
        .def_property_readonly(
            "particle_volumes",
            [](const PYDomain& self) {
                const auto values = self.get_particle_volumes();
                py::array_t<double> arr(static_cast<py::ssize_t>(values.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = values[static_cast<std::size_t>(i)];
                }

                py::object ureg = get_shared_ureg();
                return arr * ureg.attr("meter**3");
            },
            R"doc(
                Per bin particle volumes.

                Returns
                -------
                pint.Quantity[ndarray]
                    One dimensional array of per particle volumes in meter**3 computed as (4/3)π r^3.
            )doc"
        )
        .def_property_readonly(
            "total_particle_volume",
            [](const PYDomain& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_total_particle_volume()) * ureg.attr("meter**3");
            },
            R"doc(
                Total occupied particle volume.

                Returns
                -------
                pint.Quantity
                    Total occupied volume in meter**3 computed as volume_fraction * volume.
            )doc"
        )
        .def_property_readonly(
            "mean_particle_volume_number_weighted",
            [](const PYDomain& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_mean_particle_volume_number_weighted()) * ureg.attr("meter**3");
            },
            R"doc(
                Number weighted mean particle volume.

                Returns
                -------
                pint.Quantity
                    Mean per particle volume in meter**3 computed using number_fractions as weights.
            )doc"
        )
        .def_property_readonly(
            "number_of_particles_total",
            &PYDomain::get_number_of_particles_total,
            R"doc(
                Total inferred particle count.

                This is computed from:
                - total_particle_volume = volume_fraction * volume
                - mean_particle_volume_number_weighted
                and then rounded according to rounding_mode.

                Returns
                -------
                int
                    Total integer particle count.
            )doc"
        )
        .def_property_readonly(
            "number_of_particles_per_radius",
            [](const PYDomain& self) {
                const auto values = self.get_number_of_particles_per_radius();
                py::array_t<long long> arr(static_cast<py::ssize_t>(values.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = static_cast<long long>(values[static_cast<std::size_t>(i)]);
                }
                return arr;
            },
            R"doc(
                Per bin integer particle counts.

                Returns
                -------
                numpy.ndarray
                    One dimensional integer array of particle counts per radius bin.
            )doc"
        )
        .def_property_readonly(
            "particle_densities_per_radius",
            [](const PYDomain& self) {
                const auto values = self.get_particle_densities_per_radius();
                py::array_t<double> arr(static_cast<py::ssize_t>(values.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = values[static_cast<std::size_t>(i)];
                }
                py::object ureg = get_shared_ureg();
                py::object inv_m3 = py::float_(1.0) / ureg.attr("meter**3");
                return arr * inv_m3;
            },
            R"doc(
                Per bin number densities.

                Returns
                -------
                pint.Quantity[ndarray]
                    One dimensional array of number densities in 1/meter**3.
            )doc"
        )
        .def_property_readonly(
            "particle_density_total",
            [](const PYDomain& self) {
                py::object ureg = get_shared_ureg();
                py::object inv_m3 = py::float_(1.0) / ureg.attr("meter**3");
                return py::float_(self.get_particle_density_total()) * inv_m3;
            },
            R"doc(
                Total number density over all bins.

                Returns
                -------
                pint.Quantity
                    Total number density in 1/meter**3.
            )doc"
        )
        .def_property_readonly(
            "volume_fraction_per_radius",
            [](const PYDomain& self) {
                const auto values = self.get_volume_fraction_per_radius();
                py::array_t<double> arr(static_cast<py::ssize_t>(values.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = values[static_cast<std::size_t>(i)];
                }
                return arr;
            },
            R"doc(
                Per bin occupied volume fractions.

                Returns
                -------
                numpy.ndarray
                    One dimensional array of per bin volume fractions (dimensionless).
            )doc"
        )
        .def(
            "sample_radii",
            [](const PYDomain& self, long long number_of_samples, py::object seed_py) {
                std::uint64_t seed = 0;
                if (seed_py.is_none()) {
                    seed = static_cast<std::uint64_t>(std::random_device{}());
                } else {
                    seed = py::int_(seed_py);
                }

                const auto samples = self.sample_radii(number_of_samples, seed);

                py::array_t<double> arr(static_cast<py::ssize_t>(samples.size()));
                auto buf = arr.mutable_unchecked<1>();
                for (py::ssize_t i = 0; i < buf.shape(0); ++i) {
                    buf(i) = samples[static_cast<std::size_t>(i)];
                }
                py::object ureg = get_shared_ureg();
                return arr * ureg.attr("meter");
            },
            py::arg("number_of_samples"),
            py::arg("seed") = py::none(),
            R"doc(
                Sample particle radii according to number_fractions.

                Parameters
                ----------
                number_of_samples : int
                    Number of radii draws.
                seed : Optional[int]
                    RNG seed. If None, a non deterministic seed is used.

                Returns
                -------
                pint.Quantity[ndarray]
                    One dimensional Pint quantity array of sampled radii in meters.
            )doc"
        )
        .def(
            "print_bins",
            [](const PYDomain& self, int precision) { self.print_bins(precision); },
            py::arg("precision") = 6,
            R"doc(
                Print a per bin table describing the polydisperse population.

                The table is printed in base units and includes:
                - radius (m)
                - per particle volume (m^3)
                - number fraction (dimensionless)
                - volume fraction (dimensionless)
                - particle count (integer)
                - number density (1/m^3)

                Parameters
                ----------
                precision : int, default=6
                    Significant digits used for numeric formatting.

                Returns
                -------
                None
            )doc"
        )
        .def(
            "bins_table",
            &PYDomain::bins_table,
            py::arg("precision") = 6,
            R"doc(
                Return the per bin table as a string.

                Parameters
                ----------
                precision : int, default=6
                    Significant digits used for numeric formatting.

                Returns
                -------
                str
                    Formatted table in base units.
            )doc"
        )
        .def("__repr__", [](const PYDomain& self) {
            std::ostringstream stream;
            stream << "<PercusYevickDomain size=" << self.size << " m, species=" << self.radii.size() << ">";
            return stream.str();
        });
}
