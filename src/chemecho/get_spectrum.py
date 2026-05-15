import nistchempy as nist
from jcamp import jcamp_readfile
import tempfile
import os
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO
from spectrum_utils import YUNITS_PRIORITY, classify_yunits, to_transmittance, to_wavenumber


def _parse_jdx(spec) -> dict:
    """Parse a nistchempy IR spectrum object into a jcamp dict.

    Args:
        spec: a nistchempy IR spectrum object with a jdx_text attribute
    Returns:
        dict with keys 'x', 'y', 'xunits', 'yunits', etc.
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.jdx', delete=False) as f:
        f.write(spec.jdx_text)
        temp_path = f.name
    data = jcamp_readfile(temp_path)
    os.remove(temp_path)
    return data


def _get_yunits(spec) -> str:
    """Read yunits directly from jdx_text header without full parsing."""
    for line in spec.jdx_text.splitlines():
        if '##YUNITS=' in line.upper():
            return line.split('=', 1)[1].strip()
    return ''


def extract_spectrum_data(compound_or_cas):
    """Extract IR spectrum from NIST and return it as transmittance (%).

    Searches all available spectra for the compound in priority order:
        1. TRANSMITTANCE  — used directly, no conversion
        2. ABSORBANCE     — converted via T = 10^(-A) * 100
        3. MOLAR_ABSORPTIVITY — approximated as absorbance-like

    Args:
        compound_or_cas: CAS number string, or an already-fetched NistCompound object.
            Pass the compound object when you already have it to avoid a redundant
            NIST request.
    Returns:
        tuple (list, list): wavenumbers in cm⁻¹ and transmittance in %
    Raises:
        ValueError: if no spectrum in a supported format is found
    """
    if isinstance(compound_or_cas, str):
        X = nist.get_compound(compound_or_cas)
    else:
        X = compound_or_cas

    X.get_ir_spectra()

    if not X.ir_specs:
        raise ValueError(f"No IR spectra found for {X.name}")

    # Check yunits from header only, parse only the chosen spectrum
    for target in YUNITS_PRIORITY:
        for spec in X.ir_specs:
            if classify_yunits(_get_yunits(spec)) != target:
                continue
            data = _parse_jdx(spec)
            try:
                data_x = to_wavenumber(data['x'], data.get('xunits', ''))
                data_y = to_transmittance(data['y'], data.get('yunits', ''))
                return data_x, data_y
            except ValueError:
                continue

    raise ValueError(
        f"Could not read any IR spectrum for {X.name} in a supported format"
    )


def ir_graph(data: tuple, compound_name: str):
    """Plot IR spectra.

    Args:
        data (tuple): the two lists of values (wavenumbers and transmittance)
        compound_name (str): name of the compound, used for the plot title
    Returns:
        tuple (fig, ax): matplotlib figure and axes objects
    """
    df_spectrum = pd.DataFrame({"Wavenumber": data[0], "Transmittance": data[1]})

    fig, ax = plt.subplots()
    ax.plot(df_spectrum["Wavenumber"], df_spectrum["Transmittance"], color="black")
    ax.set_xlabel(r"Wavenumber [cm$^{-1}$]")
    ax.set_ylabel("Transmittance [%]")
    ax.set_title(f"IR spectrum of {compound_name}")

    ax.grid(visible=False)
    ax.set_xlim(max(data[0]), min(data[0]))
    ax.set_ylim(0, 100)

    return fig, ax


def from_df_to_csv(df_spectrum) :
    buffer_csv = StringIO()
    df_spectrum.to_csv(buffer_csv, index=False)
    return buffer_csv.getvalue()


if __name__ == "__main__":
    ir_graph(extract_spectrum_data('74-85-1'), '74-85-1')
