from PackLab.binary import interface_pint  # noqa: F401
from PackLab.units import ureg # noqa: F401

debug_mode = False  # noqa: F401


try:
    from ._version import version as __version__  # noqa: F401

except ImportError:
    __version__ = "0.0.0"

# -
