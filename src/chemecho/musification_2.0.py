import chemecho_epfl.src.chemecho.get_spectrum as spect
import numpy as np
import musicpy as mp
from musicpy.daw import *
from io import StringIO
import os


def molecular_weight_to_sound_code(compound_cas : str) -> int :

    """
    Associates an instrument to the molecular weight (heavier molecule corresponds to lower instrument)
    
    Arg:
        (str) : CAS number of compound

    Return 
        (int) : int corresponding to an instrument in General MIDI Instruments
    """
    if not isinstance(compound_cas, str):
        raise TypeError (f"Invalid type {type(compound_cas)}: CAS number must be a string")
    
    compound = spect.nist.get_compound(compound_cas)
    if compound is None:
        raise ValueError (f"Could not find a compound in the database for the specified CAS number {compound_cas}")

    molecular_weight_compound = compound.mol_weight
    if molecular_weight_compound < 50 :
        return 102 #return  "goblin" sound (light sound)
    elif molecular_weight_compound >= 50 and molecular_weight_compound < 100 :
        return 89 #return "Fantasia" sound (a little bit heavier)
    elif molecular_weight_compound >= 100 and molecular_weight_compound < 200 :
        return 92 #return "space voice" sound (intermediately light)
    elif molecular_weight_compound >= 200 and molecular_weight_compound < 300 :
        return 95 #return "Halo pad" sound (heavier)
    else :
        return 99 #return "Crystal" sound (heavy)

def peak_detection (wavenumbers, transmittances) -> list :
    
    """
    Removes the tiny noise so the final music is a little more pleasing
    
    arg (str): CAS number of compound
    return (list): list corresponding to spectrum without noise
    
    """
    
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

    The audio-spectrum frequency varies with a signal: if low transmittance:high frequency
    Saves to a MIDI file, that can be incorporated in Streamlit app (or played with VLC media player).

    Args (str), (int) : CAS number of compound, bpm of music it default is 120 bpm
    Return (str) : Filename of music generated, to be included in Streamlit
    """
    extracted_data = spect.extract_spectrum_data(cas)
    wavenumbers = extracted_data[0]
    transmittances = extracted_data[1]
    peaks = peak_detection(wavenumbers, transmittances)
    print (len(peaks))
    compound = spect.nist.get_compound(cas).name
    instru = molecular_weight_to_sound_code(cas)

    notes=[]
    target_duration_beats = 30
    interval_value = target_duration_beats / len(peaks)
    intervals = [interval_value] * len(peaks)

    for wave, trans in peaks:
        freq = 196 + (100 - trans) / 100 * (1046.50 - 196)
        note_generated = mp.freq_to_note(freq)
        note_name = note_generated.name
        octave = note_generated.num
        volume = 80  # constant volume
        note_duration = interval_value * 10 #duration longer than space btw notes that way they 
        #"mix" and the sound becomes continuous
        if freq < 200:
            note = mp.note(note_name,octave, note_duration, volume=20)
        else:
            note = mp.note(note_name,octave, note_duration, volume)
        #this is here so the "noise" is not as annoying
        note_with_effect = set_effect(note, fade(5,5))
        notes.append(note_with_effect)

  
    notes_track = mp.chord(notes=notes, interval=intervals)

    min_wave = min(wavenumbers)
    max_wave = max(wavenumbers)
    
    start_small = int(np.ceil(min_wave / 100)) * 100
    end_small = int(np.floor(max_wave / 100)) * 100
    
    print(f"interval_value: {interval_value}")
    print(f"note_duration: {note_duration}")
    print(f"durée théorique en beats: {interval_value * len(peaks)}")
    print(f"durée note en ms: {note_duration * (60000/bpm_mol)}")

    #final music composition !!!
    all_tracks = [notes_track]
    instruments = [instru]
    channels = [0]
    
    all_marks = list(range(start_small, end_small + 1, 100))
    drum_string_parts = []
    if len(all_marks) > 0: # here we check that the scale is not empty which might happen
        #(although unlikely) if a spectrum is very narrow
        for mark in all_marks:
            if mark % 500 == 0:
                drum_string_parts.append('K') # kick drum for every 500 cm-1
            else:
                drum_string_parts.append('S') # snare drum for every 100 cm-1
                # that way we don't have any 
                #overlap of the 500 cm-1 mark and the 100 cm-1 mark
    interval_drum = 30 / len(all_marks)
    drum_string = ', '.join(drum_string_parts)
    drum_track = mp.drum(drum_string, default_interval=interval_drum)

    if drum_track is not None:
        all_tracks.append(drum_track.notes)
        instruments.append(1)
        channels.append(9)


    music = mp.piece(tracks=all_tracks,
                  instruments=instruments,
                  channels=channels,
                  bpm=bpm_mol,
                  name=compound)
    filename = f"{compound}_audio_spectrum_alternative_version.mid"

    mp.write(current_chord=music, name=filename)

    return filename #returns the name of the file which can then be used in streamlit hopefully

# mini test
if __name__ == "__main__":
    result = molecular_music('57-50-1')
    print(f"File created : {result}")
    print (f"The sound chosen is: {molecular_weight_to_sound_code('57-50-1')}")