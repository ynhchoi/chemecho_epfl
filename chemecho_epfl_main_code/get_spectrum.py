import nistchempy as nist
from jcamp import jcamp_readfile
import tempfile
import os
from ir_spectra_conversion import normalize_to_transmittance
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO


def _pick_best_spec(ir_specs):
    """Pick the best spectrum: prefer TRANSMITTANCE, then ABSORBANCE, then fallback to first."""
    for preferred in ['TRANSMITTANCE', 'ABSORBANCE']:
        for spec in ir_specs:
            if f'YUNITS={preferred}' in spec.jdx_text.upper():
                return spec
    return ir_specs[0]


def extract_spectrum_data(cas: str):
    """Fetch IR spectrum data from NIST for a given CAS number.
    Returns a jcamp dict with 'x' (wavenumbers) and 'y' (transmittance 0-100%).
    """
    X = nist.get_compound(cas)
    X.get_ir_spectra()
    X_IR = _pick_best_spec(X.ir_specs)
    jdx_content = X_IR.jdx_text

    # we need a tempfile because jcamp_readfile only reads from disk
    with tempfile.NamedTemporaryFile(mode='w', suffix='.jdx', delete=False) as f:
        f.write(jdx_content)
        temp_path = f.name

    data = jcamp_readfile(temp_path)
    os.remove(temp_path)
    return normalize_to_transmittance(data)


def get_wavenumbers_and_intensities(cas: str):
    """Convenience wrapper that returns (wavenumbers, transmittance) directly."""
    data = extract_spectrum_data(cas)
    return data['x'], data['y']
    data_x = data['x']
    data_y = data['y']

    if data.get('yunits') == "Absorbance" or data.get('yunits') == "ABSORBANCE" :
        for i in range(0, len(data_y)):
            data_y[i] = 10**(-data_y[i])
    elif data.get('yunits') == "Transmittance" or data.get('yunits') == "TRANSMITTANCE" :
        pass
    else :    
        raise ValueError(f"{X.name} has an IR spectrum with y units not convertible in Transmittance")

    return data_x, data_y


def ir_graph(data : tuple, cas_number) :
    df_spectrum = pd.DataFrame({"Wavenumber":data[0], "Transmittance":data[1]})
    compound = nist.get_compound(cas_number)

    fig, ax = plt.subplots()
    ax.plot(df_spectrum["Wavenumber"], df_spectrum["Transmittance"], color="black")
    ax.set_xlabel("Wavenumber []")
    ax.set_ylabel("Transmittance [%]")
    ax.set_title(f"IR spectrum of {compound.name}")

    ax.grid(visible=False)
    ax.set_xlim(max(data[0]), min(data[0]))
    ax.set_ylim(0, 1)

    return fig.savefig("test.png")


def from_df_to_csv(df_spectrum) :
    buffer_csv = StringIO()
    return df_spectrum.to_csv(buffer_csv, index=False)

ir_graph(extract_spectrum_data('58-08-2'),'58-08-2')
