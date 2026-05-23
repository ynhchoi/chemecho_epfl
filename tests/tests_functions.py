import pytest
import sys
from pathlib import Path
import nistchempy as nist

sys.path.append(str(Path(__file__).parent.parent / "src\chemecho"))

from get_spectrum import extract_spectrum_data, ir_graph
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
    assert ir_graph(extract_spectrum_data(nist.get_compound('58-08-2')), nist.get_compound('58-08-2').name) == "<Figure size 640x480 with 1 Axes>, <Axes: title={'center': 'IR spectrum of Caffeine'}, xlabel='Wavenumber [cm$^{-1}$]', ylabel='Transmittance [%]'>"

test_extract_spectrum_data()
test_ir_graph()