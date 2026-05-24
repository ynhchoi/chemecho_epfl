"""Tests for the ``get_spectrum`` module."""

import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import nistchempy as nist
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "chemecho"))

from get_spectrum import extract_spectrum_data, ir_graph


@pytest.fixture(autouse=True)
def _nist_cooldown():
    # NIST WebBook rate-limits rapid consecutive requests; nistchempy 2.0.0
    # uses delay=0 by default so back-to-back tests can get empty responses.
    # A short cooldown between tests keeps the NIST server happy.
    yield
    time.sleep(2)


def test_extract_spectrum_data():
    # test with a CAS number that gives the first result
    result1 = extract_spectrum_data(nist.get_compound('74-85-1'))
    transmittance1 = result1[1]
    assert pytest.approx(transmittance1[0:10], abs=0.1) == [
        92.9, 92.9, 92.9, 92.9, 92.9, 92.9, 93.0, 93.0, 93.0, 93.2,
    ]
    # test with a CAS number that raises a ValueError
    with pytest.raises(
        ValueError,
        match="Could not find a spectrum for Trinitrotoluene or Trinitrotoluene has an IR spectrum with y units not convertible in Transmittance",
    ):
        extract_spectrum_data(nist.get_compound('118-96-7'))


def test_ir_graph():
    fig, ax = ir_graph(
        extract_spectrum_data(nist.get_compound('58-08-2')),
        nist.get_compound('58-08-2').name,
    )

    assert ax.get_title() == "IR spectrum of Caffeine"
    assert pytest.approx(ax.get_xlim(), abs=10) == (3966, 450)

    plt.close(fig)
