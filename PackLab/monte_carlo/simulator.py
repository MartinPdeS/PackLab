from PackLab.binary.interface_rsa import Options, Simulator
from PackLab.binary.interface_domain import Domain
from PackLab.binary.interface_radius_sampler import UniformRadiusSampler

from PackLab.monte_carlo.results import Result

class Simulator(Simulator):
    """
    Python convenience wrapper around the C++ RSA simulator.

    It creates the bound RandomSequentialAdditionSimulator, runs it, then returns an Result
    with computed statistics and plotting helpers.
    """

    def __init__(self, domain: Domain, radius_sampler: UniformRadiusSampler, options: Options) -> None:
        self.radius_sampler = radius_sampler
        self.options = options

        super().__init__(
            domain=domain,
            radius_sampler=radius_sampler,
            options=options,
        )

    def run(self) -> Result:
        """
        Run the RSA simulation.

        Returns
        -------
        Result
            The result of the simulation, including sphere positions, radii, statistics, and plotting helpers
        """
        binding = self._cpp_run()

        result = Result(
            binding=binding,
        )

        result.sphere_configuration = self.sphere_configuration

        return result

