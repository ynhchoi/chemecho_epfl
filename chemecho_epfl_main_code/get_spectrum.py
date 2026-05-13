import nistchempy as nist
from jcamp import jcamp_readfile
import tempfile
import os
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO

def extract_spectrum_data(cas:str):
    X = nist.get_compound(cas)
    X.get_ir_spectra()
    X.ir_specs
    X_IR = X.ir_specs[0]
    jdx_content = X_IR.jdx_text

    # here we need a tempfile bc the data is in a jdx that we need to go through to extract the data
    with tempfile.NamedTemporaryFile(mode='w', suffix='.jdx', delete=False) as f:
        f.write(jdx_content)
        temp_path = f.name

    data = jcamp_readfile(temp_path)
    os.remove(temp_path)

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