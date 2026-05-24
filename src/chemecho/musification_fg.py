import numpy as np
import musicpy as mp
from musicpy.daw import set_effect, fade
from rdkit import Chem
from scipy.signal import find_peaks

def noise_reduction(wavenumbers, transmittances) -> list:
    """
    Removes the tiny noise so the final music is a little more pleasing

    Any point whose transmittance is within 10 % of the spectrum's full range
    of the maximum is replaced by``max(transmittances)``, clamping the baseline 
    flat. Points below that threshold are kept as-is.

    First and last data points are removed because of possible anomalies

    Args:
        wavenumbers (list[float]): cm⁻¹ values, in any order
        transmittances (list[float]): % transmittance values aligned to
            wavenumbers

    Returns:
        list[tuple[float, float]]: one (wavenumber, transmittance) tuple
        per data point, in input order, length =len(wavenumbers) - 2
    """
    peaks = []
    threshold = max(transmittances) - 0.1 * (max(transmittances) - min(transmittances))
    for i in range(1, len(transmittances) - 1):
        if transmittances[i] < threshold:
            peaks.append((wavenumbers[i], transmittances[i]))
        else:
            peaks.append((wavenumbers[i], max(transmittances)))
    return peaks

FG_CATALOG = {
    'O-H (alcohol)':    {'smarts': '[OX2H]',           'region': (3100, 3650), 'instrument': 73, 'pitch': ('C', 6)},
    'O-H (carboxylic)': {'smarts': '[CX3](=O)[OX2H1]', 'region': (2500, 3300), 'instrument': 72, 'pitch': ('A', 5)},
    'N-H':              {'smarts': '[NX3;H1,H2]',      'region': (3250, 3500), 'instrument': 69, 'pitch': ('G', 5)},
    'C-H (sp3)':        {'smarts': '[CX4;H1,H2,H3]',   'region': (2840, 2985), 'instrument': 25, 'pitch': ('F', 5)},
    'C-H (aromatic)':   {'smarts': '[cH]',             'region': (3000, 3120), 'instrument': 26, 'pitch': ('E', 5)},
    'C#N (nitrile)':    {'smarts': 'C#N',              'region': (2200, 2260), 'instrument': 11, 'pitch': ('D', 5)},
    'C#C (alkyne)':     {'smarts': 'C#C',              'region': (2100, 2260), 'instrument': 12, 'pitch': ('D', 5)},
    'C=O (carbonyl)':   {'smarts': '[CX3]=[OX1]',      'region': (1640, 1820), 'instrument': 56, 'pitch': ('C', 5)},
    'C=C (alkene)':     {'smarts': '[CX3]=[CX3;!c]',   'region': (1600, 1680), 'instrument': 41, 'pitch': ('B', 4)},
    'aromatic ring':    {'smarts': '[a]',              'region': (1450, 1610), 'instrument': 74, 'pitch': ('A', 4)},
    'C-O':              {'smarts': '[CX4][OX2]',       'region': (1000, 1300), 'instrument': 61, 'pitch': ('G', 4)},
    'N=O (nitro)':      {'smarts': '[N+](=O)[O-]',     'region': (1300, 1600), 'instrument': 57, 'pitch': ('F', 4)},
}
# ---------------------------------------------------------------------------
# Functional group catalog.
#
#   smarts:     RDKit SMARTS pattern used on the SMILES.
#   region:     (low, high) cm-1 window in which the IR peak is expected.
#               Boundaries are slightly larger than reference values in case of
#               shifts in absorption.
#   instrument: General MIDI program number.
#   pitch:      (note name, octave) signature sound.
#
# ---------------------------------------------------------------------------

def detect_functional_groups(smiles: str) -> list:
    """
    Detect which catalog functional groups are present in a molecule.

    Args:
        smiles (str): SMILES string of the compound.

    Returns:
        list[str]: the list of FG_CATALOG entries whose SMARTS pattern matches `smiles`.
        Empty list if the SMILES cannot be parsed.
    """
    if not smiles:
        return []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    found = []
    for name, info in FG_CATALOG.items():
        patt = Chem.MolFromSmarts(info['smarts'])
        if patt is not None and mol.HasSubstructMatch(patt):
            found.append(name)
    return found


def count_carbons(smiles: str) -> int:
    """
    Count the number of carbon atoms in a molecule.

    Args:
        smiles (str): SMILES string of the compound.

    Returns:
        int: total carbon atom count. Returns 0 if ``smiles`` is empty,
        unparseable, or the molecule contains no carbon (inorganic) — in
        which case the carbon-count prelude should be skipped.    """

    if not smiles:
        return 0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)


def find_absorption_peaks(wavenumbers, transmittances, prominence_frac: float = 0.05) -> list:
    """
    Locate true local minima of transmittance (= local maxima of absorption)
    using scipy.signal.find_peaks.

    Args:
        wavenumbers (list[float])
        transmittances (list[float])
        prominence_frac (float): minimum peak prominence as a fraction of the
            transmittance range. Filters out noise; 0.05 = 5 % of full range.

    Returns:
        list[dict]: per-peak dicts with keys
            'slot'          – index into the original arrays (synced to v1 timing)
            'wavenumber'    – cm-1
            'transmittance' – % (0..100, lower means deeper absorption)
            'prominence'    – peak depth from local baseline
    """
    trans = np.asarray(transmittances, dtype=float)
    absorption = trans.max() - trans
    span = trans.max() - trans.min()
    if span <= 0:
        return []
    prominence = max(prominence_frac * span, 1e-6)
    indices, props = find_peaks(absorption, prominence=prominence)
    n = len(wavenumbers)
    out = []
    for k, i in enumerate(indices):
        if not (1 <= i <= n - 2):
            continue  # outside the music timeline
        out.append({
            'slot': int(i) - 1,
            'wavenumber': float(wavenumbers[i]),
            'transmittance': float(trans[i]),
            'prominence': float(props['prominences'][k]),
        })
    return out


def assign_peaks_to_fgs(detected_peaks: list, present_fgs: list) -> dict:
    """
    Disambiguate overlapping FG regions: every detected peak is assigned to
    the present FG whose region (a) contains the peak's wavenumber and
    (b) is the NARROWEST among the candidates. A peak in no candidate region
    is dropped (no accent fired).

    Args:
        detected_peaks (list): detected peaks
        present_fgs (list): functional groups present in compound

    Returns:
        dict[str, list[dict]]: fg_name -> peaks assigned to it.
    """
    def sort_key(fg_name: str):
        lo, hi = FG_CATALOG[fg_name]['region']
        # primary: narrower region wins; secondary: name
        return (hi - lo, fg_name)

    out = {fg: [] for fg in present_fgs}
    for peak in detected_peaks:
        w = peak['wavenumber']
        candidates = [
            fg for fg in present_fgs
            if FG_CATALOG[fg]['region'][0] <= w <= FG_CATALOG[fg]['region'][1]
        ]
        if not candidates:
            continue
        winner = min(candidates, key=sort_key)
        out[winner].append(peak)
    return {fg: hits for fg, hits in out.items() if hits}


def _scale_accent_volume(prominence: float, max_prominence: float,
                         lo: int = 60, hi: int = 127) -> int:
    """
    Map peak prominence to MIDI volume so deeper peaks sound louder.
    
    Args:
    
        prominence (float): intensity of peak
        max_prominence (float): maximum intensity
        lo (int): low volume in MIDI
        hi (int): loud vloume in MIDI
    
    Return:
    
        (int): scale
    """
    if max_prominence <= 0:
        return hi
    return int(lo + (hi - lo) * (prominence / max_prominence))


def molecular_weight_to_bpm(compound, bpm_min: int = 60, bpm_max: int = 120) -> int:
    """
    Map molecular weight to BPM via Maxwell-Boltzmann mean speed (v ∝ 1/√M).
    Lighter molecules move faster → higher BPM; heavier → lower BPM.

    The 1/√M values are linearly interpolated between two reference masses
    (M_fast=16 g/mol, M_slow=500 g/mol) and mapped onto [bpm_min, bpm_max].
    Result is clamped so the music always stays within the specified bounds
    regardless of how extreme the molecular weight is.

    Default range 60–120 BPM → music duration 15–30 s at 30 target beats.

    Args:

        compound (nist.compound.NistCompound)
        bpm_min (int): tempo of music, default is 60 bpm
        bpm_max (int): max tempo of music, default is 120 bpm

    Return:
        (int): tempo of music
    """
    import math
    M = compound.mol_weight
    M_fast, M_slow = 16.0, 500.0
    speed      = 1.0 / math.sqrt(M)
    speed_fast = 1.0 / math.sqrt(M_fast)
    speed_slow = 1.0 / math.sqrt(M_slow)
    ratio = (speed - speed_slow) / (speed_fast - speed_slow)
    ratio = max(0.0, min(1.0, ratio))
    return int(bpm_min + (bpm_max - bpm_min) * ratio)


def molecular_music_fg(extracted_data, compound, smiles: str):
    """
    IR sonification with functional-group accents, carbon count,
    tempo adapted to molecular weight, and volume varying with peak intensity.

    Args:

        extracted_data (tuple(list, list)): data extracted from spectrum
        compound (nist.compound.NistCompound) : molecule
        smiles (str): smiles of molecule
    Returns:
        tuple[str, dict]: (midi filename, legend dict).
    Raises:
        ValueError: if the spectrum has fewer than 3 points or if
            ``wavenumbers`` and ``transmittances`` have mismatched lengths.
    """
    wavenumbers = extracted_data[0]
    transmittances = extracted_data[1]
    
    if wavenumbers[0] < wavenumbers[-1]:
        wavenumbers = wavenumbers[::-1]
        transmittances = transmittances[::-1]

    if len(wavenumbers) < 3 or len(wavenumbers) != len(transmittances):
        raise ValueError(
            "Spectrum must have at least 3 points and matching wavenumber/"
            f"transmittance lengths (got {len(wavenumbers)} / {len(transmittances)})."
        )
    peaks = noise_reduction(wavenumbers, transmittances)
    compound_name = compound.name
    instru_main = 89
    bpm_mol = molecular_weight_to_bpm(compound)
    n_slots = len(peaks)

    #carbons prelude
    n_carbons = count_carbons(smiles)
    intro_hit_interval = 0.3   
    intro_hit_duration = 0.2  
    intro_pause        = 0.2  
    intro_duration_beats = (
        n_carbons * intro_hit_interval + intro_pause if n_carbons > 0 else 0
    )

    #main melody track
    target_duration_beats = 30
    interval_value = target_duration_beats / n_slots
    intervals = [interval_value] * n_slots
    note_duration = interval_value * 8

    notes = []
    for wave, trans in peaks:
        freq = 196 + (100 - trans) / 100 * (1046.50 - 196)
        gen = mp.freq_to_note(freq)
        vol = 20 if freq < 200 else 80
        n = mp.note(gen.name, gen.num, note_duration, volume=vol)
        notes.append(set_effect(n, fade(5, 5)))
    notes_track = mp.chord(notes=notes, interval=intervals)

    #peaks detection
    present_fgs = detect_functional_groups(smiles)
    detected = find_absorption_peaks(wavenumbers, transmittances)
    confirmed = assign_peaks_to_fgs(detected, present_fgs)

    max_prom = max(
        (p['prominence'] for hits in confirmed.values() for p in hits),
        default=0.0,
    )

    #tracks assembly
    all_tracks = [notes_track]
    instruments = [instru_main]
    channels = [0]
    start_times = [intro_duration_beats]

    if n_carbons > 0:
        intro_notes = [mp.note('C', 3, intro_hit_duration, volume=110)
                       for _ in range(n_carbons)]
        intro_intervals = [intro_hit_interval] * n_carbons
        intro_track = mp.chord(notes=intro_notes, interval=intro_intervals)
        all_tracks.append(intro_track)
        instruments.append(119)  # 119 = Synth Drum
        channels.append(1)
        start_times.append(0)
        next_channel = 2         
    else:
        next_channel = 1
    

    accent_duration = interval_value * 9  
    for fg, hits in confirmed.items():
        info = FG_CATALOG[fg]
        pname, poct = info['pitch']
        hits_sorted = sorted(hits, key=lambda p: p['slot'])
        accent_notes = []
        accent_intervals = []
        for i, peak in enumerate(hits_sorted):
            vol = _scale_accent_volume(peak['prominence'], max_prom)
            an = mp.note(pname, poct, accent_duration, volume=vol)
            accent_notes.append(set_effect(an, fade(2, 10)))
            if i < len(hits_sorted) - 1:
                gap = (hits_sorted[i + 1]['slot'] - peak['slot']) * interval_value
                accent_intervals.append(gap)
            else:
                accent_intervals.append(accent_duration)
        accent_track = mp.chord(notes=accent_notes, interval=accent_intervals)
        all_tracks.append(accent_track)
        instruments.append(info['instrument'])
        start_times.append(intro_duration_beats + hits_sorted[0]['slot'] * interval_value)
        if next_channel == 9:
            next_channel = 10
        channels.append(next_channel)
        next_channel += 1

    #drum markers every 100 / 500 cm-1
    min_wave = min(wavenumbers)
    max_wave = max(wavenumbers)
    start_small = int(np.ceil(min_wave / 100)) * 100
    end_small = int(np.floor(max_wave / 100)) * 100
    all_marks = list(range(start_small, end_small + 1, 100))
    if all_marks:
        drum_parts = ['K' if m % 500 == 0 else 'S' for m in all_marks]
        interval_drum = 30 / len(all_marks)
        drum_track = mp.drum(', '.join(drum_parts), default_interval=interval_drum)
        all_tracks.append(drum_track.notes)
        instruments.append(1)
        channels.append(9)
        start_times.append(intro_duration_beats)

    music = mp.piece(tracks=all_tracks,
                     instruments=instruments,
                     channels=channels,
                     start_times=start_times,
                     bpm=bpm_mol,
                     name=compound_name)
    filename = f"{compound_name}_audio_spectrum_fg.mid"
    mp.write(current_chord=music, name=filename)

    legend = {
        'fgs': {
            fg: {
                'instrument_code': FG_CATALOG[fg]['instrument'],
                'region': FG_CATALOG[fg]['region'],
                'pitch': FG_CATALOG[fg]['pitch'],
                'n_peaks': len(confirmed[fg]),
            }
            for fg in confirmed
        },
        'carbon_count': n_carbons,
        'bpm': bpm_mol,
    }
    return filename, legend