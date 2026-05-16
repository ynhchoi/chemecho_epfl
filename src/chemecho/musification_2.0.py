import nistchempy as nist
from midiutil import MIDIFile
import math
from scipy.signal import find_peaks
from get_spectrum import extract_spectrum_data


def mw_to_bpm(mol_weight: float) -> int:
    """Molecular weight (Da) → tempo (BPM).

    Based on kinetic theory: RMS molecular velocity scales as 1/√M,
    so lighter molecules get a faster tempo.
    Range: 18 Da (water) → 170 BPM, 600 Da (large drug) → 60 BPM.

    Args:
        mol_weight: molecular weight in Da
    Returns:
        tempo in beats per minute (60–170)
    """
    ratio = (math.sqrt(mol_weight) - math.sqrt(18)) / (math.sqrt(600) - math.sqrt(18))
    return int(max(60, min(170, round(170 - ratio * 110))))


def molecular_weight_to_sound_code(compound) -> int:
    """Associates an instrument to the molecular weight.

    Heavier molecule corresponds to a lower-pitched instrument,
    following the analogy that heavier objects vibrate more slowly.

    Args:
        compound: NistCompound object (already fetched, avoids redundant NIST call)
    Returns:
        General MIDI program number (0-indexed)
    """
    mw = compound.mol_weight
    if mw < 50:
        return 40   # violin
    elif mw < 100:
        return 71   # clarinet
    elif mw < 200:
        return 42   # cello
    elif mw < 300:
        return 66   # tenor sax
    else:
        return 43   # contrabass


def peak_detection(wavenumbers, transmittances,
                   min_prominence: float = 15.0, min_width_pts: int = 3):
    """Detect IR absorption peaks as local minima in transmittance.

    Uses scipy find_peaks on the inverted transmittance curve to locate
    absorption dips. Returns wavenumber, prominence (absorption depth),
    and FWHM half-width for each peak.

    Args:
        wavenumbers: list of wavenumbers in cm⁻¹
        transmittances: list of transmittance values in %
        min_prominence: minimum absorption depth (%) to qualify as a peak.
                        Filters out baseline noise.
        min_width_pts: minimum peak width in data points.
                        Filters out single-sample noise spikes.
    Returns:
        list of (wavenumber, prominence, width_wn) tuples,
        sorted by wavenumber descending (high to low)
    """
    neg_t = [-v for v in transmittances]
    indices, props = find_peaks(neg_t, prominence=min_prominence, width=min_width_pts)

    # Convert width from data points to wavenumber units
    avg_spacing = abs(wavenumbers[-1] - wavenumbers[0]) / max(len(wavenumbers) - 1, 1)

    peaks = []
    for idx, prom, width_pts in zip(indices, props['prominences'], props['widths']):
        width_wn = float(width_pts) * avg_spacing
        peaks.append((wavenumbers[idx], float(prom), float(width_wn)))

    peaks.sort(key=lambda p: p[0], reverse=True)
    return peaks


def wavenumber_to_midi(wavenumber: float) -> int:
    """Map wavenumber (cm⁻¹) to MIDI pitch using linear scaling.

    Higher wavenumber = faster molecular vibration = higher pitch.
    Range: 400 cm⁻¹ → MIDI 40 (E2), 4000 cm⁻¹ → MIDI 84 (C6).

    Args:
        wavenumber: IR wavenumber in cm⁻¹
    Returns:
        MIDI pitch number (40–84)
    """
    ratio = (wavenumber - 400) / (4000 - 400)
    return int(round(max(40, min(84, 40 + ratio * 44))))


def molecular_music(compound_or_cas, output_path: str, data=None,
                    min_prominence: float = 15.0, min_width_pts: int = 3,
                    drum_track: bool = True):
    """Translate an IR spectrum into a MIDI file.

    Each detected absorption peak becomes one note, with four independent
    dimensions of spectral information encoded:

        wavenumber gap → note interval  (peaks close in wn play rapidly;
                                         peaks far apart have a longer pause)
        wavenumber     → pitch          (linear: 400 cm⁻¹ = E2, 4000 cm⁻¹ = C6)
        prominence     → velocity       (stronger absorption = louder note)
        FWHM width     → duration       (broader peak = longer, sustained note)

    Total melody duration is normalised to 30 beats so all molecules produce
    a similarly-lengthed piece while preserving relative spectral spacing.

    An optional drum track marks every 100 cm⁻¹ (hi-hat) and every
    500 cm⁻¹ (kick drum) as a spectral ruler across the full duration.

    Tempo and instrument are derived from molecular weight via kinetic theory.

    Args:
        compound_or_cas: CAS string or already-fetched NistCompound object
        output_path: file path to write the .mid file
        data: optional pre-fetched (wavenumbers, transmittances) tuple.
              Pass this to avoid fetching from NIST a second time.
        min_prominence: minimum absorption depth (%) to count as a peak
        min_width_pts: minimum peak width in data points
        drum_track: if True, add a percussion track as a wavenumber ruler

    Example:
        >>> molecular_music('67-64-1', 'acetone.mid')
    """
    if isinstance(compound_or_cas, str):
        compound = nist.get_compound(compound_or_cas)
    else:
        compound = compound_or_cas

    if data is not None:
        wavenumbers, transmittances = data
    else:
        wavenumbers, transmittances = extract_spectrum_data(compound)

    peaks = peak_detection(wavenumbers, transmittances, min_prominence, min_width_pts)
    if not peaks:
        raise ValueError(
            f"No peaks detected for {compound.name} with current thresholds. "
            "Try lowering min_prominence."
        )

    bpm        = mw_to_bpm(compound.mol_weight)
    instrument = molecular_weight_to_sound_code(compound)

    target_beats = 30.0

    # Note intervals: proportional to wavenumber gaps between consecutive peaks.
    # Peaks close in wavenumber space → short interval (rapid succession).
    # Peaks far apart → longer interval (silence between them).
    # Intervals are normalised so their sum equals target_beats.
    wns = [p[0] for p in peaks]
    if len(peaks) > 1:
        gaps = [wns[i] - wns[i + 1] for i in range(len(wns) - 1)]
        min_gap, max_gap = min(gaps), max(gaps)
        gap_range = max(max_gap - min_gap, 1.0)
        raw_intervals = [0.5 + (g - min_gap) / gap_range * 3.5 for g in gaps]
        raw_intervals.append(raw_intervals[-1])  # last note reuses previous interval
        scale = target_beats / sum(raw_intervals)
        intervals = [iv * scale for iv in raw_intervals]
    else:
        intervals = [target_beats]

    # Note duration: proportional to FWHM peak width (0.5–2.0 beats).
    # Narrower peaks → short note; broader peaks → sustained note.
    # Independent of interval so spacing and sustain encode separate information.
    all_widths  = [w for _, _, w in peaks]
    min_w, max_w = min(all_widths), max(all_widths)
    width_range  = max(max_w - min_w, 1.0)

    n_tracks = 2 if drum_track else 1
    midi = MIDIFile(n_tracks)
    midi.addTempo(0, 0, bpm)
    midi.addProgramChange(0, 0, 0, instrument)

    current_beat = 0.0
    for i, (wavenumber, prominence, width_wn) in enumerate(peaks):
        pitch       = wavenumber_to_midi(wavenumber)
        velocity    = int(max(45, min(110, 45 + (prominence / 100) * 65)))
        width_ratio = (width_wn - min_w) / width_range
        duration    = 0.5 + width_ratio * 1.5      # 0.5–2.0 beats from FWHM
        midi.addNote(0, 0, pitch, current_beat, duration, velocity)
        current_beat += intervals[i]

    if drum_track:
        midi.addTempo(1, 0, bpm)
        start_mark = int(math.ceil(min(wavenumbers) / 100)) * 100
        end_mark   = int(math.floor(max(wavenumbers) / 100)) * 100
        all_marks  = list(range(start_mark, end_mark + 1, 100))

        if all_marks:
            drum_interval = target_beats / len(all_marks)
            for j, mark in enumerate(all_marks):
                # Kick drum (MIDI 36) every 500 cm⁻¹, hi-hat (MIDI 42) every 100 cm⁻¹
                drum_pitch = 36 if mark % 500 == 0 else 42
                midi.addNote(1, 9, drum_pitch, j * drum_interval, drum_interval * 0.5, 80)

    with open(output_path, "wb") as f:
        midi.writeFile(f)
