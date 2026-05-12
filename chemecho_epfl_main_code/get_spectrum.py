import nistchempy as nist
from jcamp import jcamp_readfile
import tempfile
import os

def extract_spectrum_data(cas:str):
    X = nist.get_compound(cas)
    X.get_ir_spectra()
    X_IR = X.ir_specs[0]
    jdx_content = X_IR.jdx_text

    # here we need a tempfile bc the data is in a jdx that we need to go through to extract the data
    with tempfile.NamedTemporaryFile(mode='w', suffix='.jdx', delete=False) as f:
        f.write(jdx_content)
        temp_path = f.name

    data = jcamp_readfile(temp_path)
    os.remove(temp_path)

    # extract data
    wavenumbers = data['x']
    intensities = data['y']
    return data

print(extract_spectrum_data('74-85-1').get('yunits'))
print("3")