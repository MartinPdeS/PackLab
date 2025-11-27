"""
PackLab test.
"""

import pytest
import PackLab


def test_dummy():
    assert 1 == 1  # Simple sanity check


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
