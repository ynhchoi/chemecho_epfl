"""
musification Layer 2 — functional-group-aware sonification (approach C, hybrid).

Built on top of Layer 1 (musification.py). The main melody track is the same as L1
(pitch follows the transmittance contour). On top of it, for each functional
group that is BOTH present in the SMILES AND confirmed by a real absorption
peak in its canonical IR region, an accent track plays a signature instrument
on a signature pitch at the matching peak position.

Listener mapping:
    contour (main track)  -> overall IR spectrum shape
    accent instrument     -> which functional group
    accent timing         -> where the corresponding IR peak sits
    accent loudness       -> peak depth / prominence
"""

import numpy as np
import musicpy as mp
from musicpy.daw import set_effect, fade
from rdkit import Chem
from scipy.signal import find_peaks

from musification import molecular_weight_to_sound_code, peak_detection


# ---------------------------------------------------------------------------
# Functional group catalog.
#
#   smarts:     RDKit SMARTS pattern used on the SMILES.
#   region:     (low, high) cm-1 window in which the IR peak is expected.
#               Boundaries are slightly widened over textbook ranges to absorb
#               conjugation, H-bonding, and ring-strain shifts (e.g. amide
#               C=O can drop to ~1650).
#   instrument: General MIDI program number.
#   pitch:      (note name, octave) signature pitch.
#
# When a detected peak falls inside multiple regions belonging to FGs that are
# ALL present in the molecule, the FG with the NARROWEST region wins
# (most-specific assignment) — see assign_peaks_to_fgs. This prevents, for
# instance, a 2900 cm-1 C-H sp3 peak from being mis-claimed by the broad
# O-H (carboxylic) window (2500-3300).
# ---------------------------------------------------------------------------

FG_CATALOG = {
    'O-H (alcohol)':    {'smarts': '[OX2H]',           'region': (3100, 3650), 'instrument': 73, 'pitch': ('C', 6)},
    'O-H (carboxylic)': {'smarts': '[CX3](=O)[OX2H1]', 'region': (2500, 3300), 'instrument': 71, 'pitch': ('A', 5)},
    'N-H':              {'smarts': '[NX3;H1,H2]',      'region': (3250, 3500), 'instrument': 68, 'pitch': ('G', 5)},
    'C-H (sp3)':        {'smarts': '[CX4;H1,H2,H3]',   'region': (2840, 2985), 'instrument': 25, 'pitch': ('F', 5)},
    'C-H (aromatic)':   {'smarts': '[cH]',             'region': (3000, 3120), 'instrument': 26, 'pitch': ('E', 5)},
    'C#N (nitrile)':    {'smarts': 'C#N',              'region': (2200, 2260), 'instrument': 11, 'pitch': ('D', 5)},
    'C#C (alkyne)':     {'smarts': 'C#C',              'region': (2100, 2260), 'instrument': 12, 'pitch': ('D', 5)},
    'C=O (carbonyl)':   {'smarts': '[CX3]=[OX1]',      'region': (1640, 1820), 'instrument': 56, 'pitch': ('C', 5)},
    'C=C (alkene)':     {'smarts': '[CX3]=[CX3;!c]',   'region': (1600, 1680), 'instrument': 40, 'pitch': ('B', 4)},
    'aromatic ring':    {'smarts': '[a]',              'region': (1450, 1610), 'instrument': 74, 'pitch': ('A', 4)},
    'C-O':              {'smarts': '[CX4][OX2]',       'region': (1000, 1300), 'instrument': 60, 'pitch': ('G', 4)},
    'N=O (nitro)':      {'smarts': '[NX3](=O)=O',      'region': (1300, 1600), 'instrument': 57, 'pitch': ('F', 4)},
}


def detect_functional_groups(smiles: str) -> list:
    """
    Return the list of FG_CATALOG entries whose SMARTS pattern matches `smiles`.
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
    Total carbon atom count in the molecule. Returns 0 if SMILES is empty,
    invalid, or the molecule contains no carbon (inorganic) — in which case
    the carbon-count intro should be skipped.
    """
    if not smiles:
        return 0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)


def find_absorption_peaks(wavenumbers, transmittances, prominence_frac: float = 0.05) -> list:
    """
    Locate true local minima of transmittance (= local maxima of absorption)
    using scipy.signal.find_peaks. Unlike v1's threshold-based peak_detection
    (which marks every slot below a threshold and therefore fires accents on
    every slot in a broad band), this returns ONE entry per actual peak.

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
    # Layer 1's peak_detection skips the first and last samples (range(1, len-1)),
    # so the music timeline uses indices 1..len-2 of the original arrays.
    # We shift slot by -1 here so it directly indexes the music timeline,
    # which keeps accent timing in sync with the main melody track.
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

    Returns:
        dict[str, list[dict]]: fg_name -> peaks assigned to it.
    """
    def sort_key(fg_name: str):
        lo, hi = FG_CATALOG[fg_name]['region']
        # primary: narrower region wins; secondary: name (deterministic tie-break)
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
    """Map peak prominence to MIDI velocity so deeper peaks sound louder."""
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
    Hybrid IR sonification with functional-group accents.

    Improvements over the first version of this file:
      - True local-minimum peak detection (scipy.signal.find_peaks) replaces
        threshold filtering. Each absorption band now fires ONE accent at its
        deepest point, not a dense burst across the whole band.
      - Narrowness-based disambiguation when FG regions overlap, so a peak
        is never claimed by two functional groups simultaneously.
      - FG regions slightly widened to catch conjugated / H-bonded shifts
        (e.g. amide C=O at ~1660).
      - "aromatic ring" SMARTS uses '[a]' so heteroaromatics (pyridine,
        furan, imidazole, ...) also trigger the aromatic accent.
      - Accent volume scales with peak prominence (deeper peak = louder),
        preserving the dynamic information lost in the constant-volume design.

    Returns:
        tuple[str, dict]: (midi filename, legend dict).
    """
    wavenumbers = extracted_data[0]
    transmittances = extracted_data[1]
    if len(wavenumbers) < 3 or len(wavenumbers) != len(transmittances):
        raise ValueError(
            "Spectrum must have at least 3 points and matching wavenumber/"
            f"transmittance lengths (got {len(wavenumbers)} / {len(transmittances)})."
        )
    peaks = peak_detection(wavenumbers, transmittances)
    compound_name = compound.name
    instru_main = 42  # Cello
    bpm_mol = molecular_weight_to_bpm(compound)
    n_slots = len(peaks)

    # --- carbon-count prelude ---
    # For organic molecules, play one Synth Drum hit per carbon at the start,
    # quickly, followed by a short pause before the spectrum begins. Same
    # pitch (C3) for every hit so the listener counts the number, not a
    # melody. Inorganic molecules (no carbon) skip the prelude entirely.
    #
    # Unit note: musicpy's `duration` and `interval` parameters are measured
    # in whole notes (1.0 musicpy unit = 4 MIDI beats). The values below are
    n_carbons = count_carbons(smiles)
    intro_hit_interval = 0.3   
    intro_hit_duration = 0.2  
    intro_pause        = 0.2  
    intro_duration_beats = (
        n_carbons * intro_hit_interval + intro_pause if n_carbons > 0 else 0
    )

    # --- 1) main melody track (v1 behavior preserved verbatim) ---
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

    # --- 2) detect real peaks; SMARTS-presence × peak-confirmation × disambiguation ---
    present_fgs = detect_functional_groups(smiles)
    detected = find_absorption_peaks(wavenumbers, transmittances)
    confirmed = assign_peaks_to_fgs(detected, present_fgs)

    max_prom = max(
        (p['prominence'] for hits in confirmed.values() for p in hits),
        default=0.0,
    )

    # --- 3) tracks assembly: main melody + (intro) + accents + drum ---
    # All non-intro tracks are shifted by intro_duration_beats so the prelude
    # plays first, then the spectrum starts. Channel 1 is reserved for the
    # Tubular Bells intro when present.
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
        instruments.append(119)  # GM #119 = Synth Drum (short, snappy)
        channels.append(1)
        start_times.append(0)
        next_channel = 2         # accents start at channel 2 when intro present
    else:
        next_channel = 1
    

    accent_duration = interval_value * 8  
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

    # --- 4) drum markers every 100 / 500 cm-1 (unchanged from v1) ---
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


# mini test
if __name__ == "__main__":
    import nistchempy as nist
    from get_spectrum import extract_spectrum_data

    # acetone: CAS 67-64-1, SMILES "CC(=O)C"
    cas = '67-64-1'
    smiles = 'CC(=O)C'
    compound = nist.get_compound(cas)
    fname, legend = molecular_music_fg(extract_spectrum_data(compound), compound, smiles)
    print(f"File created: {fname}")
    print(f"Prelude: {legend['carbon_count']} Synth Drum hits (carbon count)")
    print(f"BPM: {legend['bpm']}  (MW-derived, Maxwell-Boltzmann)")
    print("Functional-group legend:")
    for fg, info in legend['fgs'].items():
        print(f"  {fg}: instrument={info['instrument_code']}, "
              f"region={info['region']} cm-1, pitch={info['pitch']}, "
              f"n_peaks={info['n_peaks']}")
