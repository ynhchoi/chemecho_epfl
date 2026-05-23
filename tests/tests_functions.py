import pytest
import sys
from pathlib import Path
import nistchempy as nist
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).parent.parent / "src\chemecho"))

from get_spectrum import extract_spectrum_data, ir_graph
from streamlit_app_v2 import name_to_cas, get_smiles, draw_molecule
#from musification import 

def test_extract_spectrum_data():
    #test with a CAS number that give the first result
    result1 = extract_spectrum_data(nist.get_compound('74-85-1'))
    transmittance1 = result1[1]
    assert pytest.approx(transmittance1[0:10],abs=0.1) == [92.9, 92.9, 92.9, 92.9, 92.9, 92.9, 93.0, 93.0, 93.0, 93.2]
    #test with a CAS number that raise a ValueError
    with pytest.raises(ValueError, match="Could not find a spectrum for Trinitrotoluene or Trinitrotoluene has an IR spectrum with y units not convertible in Transmittance"):
        extract_spectrum_data(nist.get_compound('118-96-7'))
    #test with a CAS number that use the recursive call of the function
    #result2 = extract_spectrum_data('')[1]
    #assert pytest.approx(result2[0,10],abs=0.01) == []

def test_ir_graph():
    fig, ax = ir_graph(extract_spectrum_data(nist.get_compound('58-08-2')), nist.get_compound('58-08-2').name)

    assert ax.get_title() == "IR spectrum of Caffeine"
    assert pytest.approx(ax.get_xlimit(),abs=10) == 4000, 450

    plt.close(fig)

def test_name_to_cas():
    #tests with two different molecules
    assert name_to_cas("acetone") == "67-64-1"
    assert name_to_cas("benzene") == "71-43-2"
    assert name_to_cas("cholesterol") == "57-88-5"
    #test with a molecule name that is not found on NIST database
    with pytest.raises(ValueError, match="No compound with IR data found for 'quinic acid' on NIST.")
        name_to_cas("quinic acid")

def test_get_smiles():
    #tests with two different molecules
    assert get_smiles("ethylene") == "C=C"
    assert get_smiles("ethanol") == "CCO"
    #test with an empty sequence, should return "None"
    assert get_smiles("") == None

def test_draw_molecule():
    #test the size of the image
    assert draw_molecule().size == 400, 300
    #test with an empty sequence, should return "None"
    assert draw_molecule("") == None


test_extract_spectrum_data()
test_ir_graph()