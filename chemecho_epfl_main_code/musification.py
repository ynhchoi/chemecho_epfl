import get_spectrum as spect
import numpy as np
import musicpy as mp
from musicpy.daw import *
from io import StringIO
import os

def wavenumber_to_frequency(wavenumber : float) -> float:

    """
    Converts wavenumber to frequency for sound in music.
    
    arg (float): wavenumber
    return (float) : frequency
    """

    frequency = 312.5*(wavenumber*0.001)**2
    return frequency

def molecular_weight_to_sound_code(compound_cas : str) -> int :

    """
    Associates an instrument to the molecular weight (heavier molecule corresponds to lower instrument)
    
    arg (str) : CAS number of compound
    return (int) : int corresponding to an instrument in General MIDI Instruments
    """

    compound = spect.nist.get_compound(compound_cas)
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
        return 59 #return tuba

def peak_detection (cas) : #je pense que je vais devoir modif les arguments
    
    """
    Removes the tiny noise so the final music is a little more pleasing
    
    arg (str): CAS number of compound
    return (list): list corresponding to spectrum without noise
    
    """
    
    wavenumbers, transmittances = spect.extract_spectrum_data(cas)
    peaks =[]
    threshold = max(transmittances)-0.1*(max(transmittances)-min(transmittances))
    for i in range(1, len(transmittances) - 1):
        if transmittances[i] < threshold:
            peaks.append((wavenumbers[i], transmittances[i]))
        else:
            peaks.append((wavenumbers[i], max(transmittances)))
    return peaks
    

def molecular_music (cas, bpm_mol=120):

    """
    Translates IR spectrum to music.

    The audio-spectrum frequency decreases with wavenumber, and the volume of a peak corresponds to its intensity.
    Saves to a MIDI file, that can be incorporated in Streamlit app (or played with VLC media player).

    args (str), (int) : CAS number of compound, bpm of music it default is 120 bpm
    return (str) : Filename of music generated, to be included in Streamlit
    """
    extracted_data = spect.extract_spectrum_data(cas)
    wavenumbers = extracted_data[0]
    transmittances = extracted_data[1]
    peaks = peak_detection(cas)
    compound = spect.nist.get_compound(cas).name
    instru = molecular_weight_to_sound_code(cas)

    notes=[]
    duration = 60/bpm_mol
    for wave, trans in peaks:
        freq = wavenumber_to_frequency(wave) #convert wavenumbers to frequencies
        note_generated = mp.freq_to_note(freq) # convert frequency to note
        note_name = note_generated.name
        octave = note_generated.num
        volume = max(5, min(127, int(((100 - trans)*1.27)))) #volume of note varies with tansmittance
        note_duration = duration * 1.2 #duration longer than space btw notes that way they 
        #"mix" and the sound becomes continuous
        note = mp.note(note_name,octave, note_duration, volume)
        note_with_effect = set_effect(note, fade(500,500))
        notes.append(note_with_effect)
    notes_track = mp.track(content=notes, instrument=instru, bpm=bpm_mol)
    
    
    #big scale for a sound every 500 cm^-1
    big_scale = []
    min_wave = min(wavenumbers)
    max_wave = max(wavenumbers)
    start_scale = int(np.ceil(min_wave / 500)) * 500 # finds first multiple of 500 that way we 
    # don't have a random scale
    end_scale = int(np.floor(max_wave / 500)) * 500 # same thing for the end
    big_scale = list(range(start_scale, end_scale + 1, 500))

    #caclulation of the total time of the music bc it is needed for the drums
    total_duration_seconds = len(peaks) * 1.2 * (60 / bpm_mol)
    measure_duration = (60 / bpm_mol) * 4
    total_measures = total_duration_seconds / measure_duration
    
    # here we check that the scale is not empty which might happen (although unlikely) 
    # if a spectrum is very narrow
    if len(big_scale) > 0:
        interval = total_measures / len(big_scale)
        big_scale_string = ', '.join(['K'] * len(big_scale))
        big_scale_track = mp.drum(big_scale_string, default_interval=interval)
    else:
        big_scale_track = None
    

    # small scale for a sound every 100 cm^-1
    small_scale = []
    start_small = int(np.ceil(min_wave / 100)) * 100
    end_small = int(np.floor(max_wave / 100)) * 100
    all_small = list(range(start_small, end_small + 1, 100))
    small_scale = [s for s in all_small if s % 500 != 0]  # removes the multiples of 500 
    #so we don't have any overlap
    

    if len(small_scale) > 0:
        interval_small = total_measures / len(small_scale)
        small_scale_string = ', '.join(['S'] * len(small_scale))
        small_scale_track = mp.drum(small_scale_string, default_interval=interval_small)
    else:
        small_scale_track = None


    #final music composition !!!
    all_tracks = [notes_track]
    instruments = [instru]
    channels = [0]
    
    if big_scale_track is not None:
        all_tracks.append(big_scale_track.notes)
        instruments.append(0)
        channels.append(9)
    
    if small_scale_track is not None:
        all_tracks.append(small_scale_track.notes)
        instruments.append(0)
        channels.append(9)

    music = mp.piece(tracks=all_tracks,
                  instruments=instruments,
                  channels=channels,
                  bpm=bpm_mol,
                  name=compound)
    filename = f"{compound}_audio_spectrum.mid"

    mp.write(current_chord=music, name=filename)

    return filename #returns the name of the file which can then be used in streamlit 
    #with st.audio(vriable contraining the result of the function, format="audio/midi")

if __name__ == "__main__":
    # mini test
    result = molecular_music("58-08-2")
    print(f"Fichier créé : {result}")