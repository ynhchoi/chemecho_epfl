import get_spectrum as spect
import numpy as np
import musicpy as mp
from io import StringIO

def wavenumber_to_frequency(wavenumber : float) -> float:
    frequency = 312.5*(wavenumber*0.001)**2
    return frequency

def molecular_weight_to_sound_code(compound_cas : str) -> int :
    compound = nist.get_compound(compound_cas)
    molecular_weight_compound = compound.mol_weight
    if molecular_weight_compound < 50 :
        return 41 #return violin sound
    elif molecular_weight_compound >= 50 and molecular_weight_compound < 100 :
        return 72 #return clarinet sound
    elif molecular_weight_compound >= 100 and molecular_weight_compound < 200 :
        return 43 #return cello sound
    elif molecular_weight_compound >= 200 and molecular_weight_compound < 300 :
        return 67 #return tenor sax sound
    else :
        return 44 #return contrabass

def peak_detection (frequencies, transmittances, threshold=0.95) :
    
    peaks =[]
    for i in range(1, len(transmittances) - 1):
        if transmittances[i] > threshold:
            peaks.append((frequencies[i], transmittances[i]))
    return peaks
    

def molecular_music (cas, bpm=120, threshold=0.95):

print("4")
