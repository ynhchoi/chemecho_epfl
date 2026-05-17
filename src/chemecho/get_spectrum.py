import nistchempy as nist
from jcamp import jcamp_readfile
import tempfile
import os
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO


def extract_spectrum_data(compound, index_spectrum: int = 0):
    """extract .jdx files of IR spectrum from NIST documentation

    Args:
        cas (str): CAS number of the compound of interest
        index_spectrum (int): the index of the spectra wanted in the list of spectrum of the compound of interest

    Return:
        tuple (list, list): two lists of values, the wavenumbers in 1/cm and the transmittance in %
    """

    #if not isinstance(cas, str):
    #    raise TypeError (f"Invalid type {type(cas)}: CAS number must be a string")

    
    if compound is None:
        raise ValueError (f"Could not find a compound in the database for the specified CAS number {compound.cas_rn}")
  
    compound.get_ir_spectra()
    compound.ir_specs

    if len(compound.ir_specs) - index_spectrum == 0:
        raise ValueError(f"Could not find a spectrum for {compound.name} or {compound.name} has an IR spectrum with y units not convertible in Transmittance")
    else: 
        compound_IR = compound.ir_specs[index_spectrum]
        jdx_content = compound_IR.jdx_text

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
            max_transmittance = max(data_y)
            for i in range(0, len(data_y)):
                data_y[i] = (data_y[i]/max_transmittance)*100
        elif data.get('yunits') == "Transmittance" or data.get('yunits') == "TRANSMITTANCE" :
            max_transmittance = max(data_y)
            for i in range(0, len(data_y)):
                data_y[i] = (data_y[i]/max_transmittance)*100
        else :    
            return extract_spectrum_data(compound,index_spectrum+1)

        if data.get('xunits') == "MICROMETERS" or data.get('xunits') == "Micrometers":
            for i in range(0, len(data_x)):
                data_x[i] = 10000/data_x[i]
        elif data.get('xunits') == "1/CM" or data.get('xunits') == "1/cm" or data.get('xunits') == "cm-1" or data.get('xunits') == "cm -1" :
            pass
        else :
            return extract_spectrum_data(compound,index_spectrum+1)

        return data_x, data_y


def ir_graph(data : tuple, compound_name : str) :
    """Plot IR spectra in .svg format

    Args:
        data (tuple): the two lists of values (wavenumbers and transmittance)
        compound_name (str): name of the compound
    Returns:
        .svg: IR spectra
    """
    df_spectrum = pd.DataFrame({"Wavenumber":data[0], "Transmittance":data[1]})

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
    return df_spectrum.to_csv(buffer_csv, index=False)

#ir_graph(extract_spectrum_data('50-78-2'),'50-78-2')
"""if __name__ == "__main__":
    result_graph = ir_graph(extract_spectrum_data('50-78-2'),'50-78-2')
    print(f"IR spectrum plotted for 50-78-2")"""
