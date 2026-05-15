import nistchempy as nist
from midiutil import MIDIFile
from get_spectrum import extract_spectrum_data


def mw_to_bpm(mol_weight: float) -> int:
    """Molecular weight (Da) → tempo (BPM).

    Lighter molecules move faster (kinetic theory analogy).
    Range: 18 Da (water) → 170 BPM, 600 Da (large drug) → 60 BPM.

    Args:
        mol_weight: molecular weight in Da
    Returns:
        tempo in beats per minute (60–170)
    """
    bpm = 170 - ((mol_weight - 18) / (600 - 18)) * 110
    return int(max(60, min(170, bpm)))


def mw_to_instrument(mol_weight: float) -> int:
    """Molecular weight → General MIDI instrument number.

    Lighter molecules → higher-pitched instruments (violin),
    heavier molecules → lower-pitched instruments (contrabass).

    Args:
        mol_weight: molecular weight in Da
    Returns:
        General MIDI program number (0-indexed)
    """
    if mol_weight < 50:
        return 40   # violin
    elif mol_weight < 100:
        return 71   # clarinet
    elif mol_weight < 200:
        return 42   # cello
    elif mol_weight < 300:
        return 66   # tenor sax
    else:
        return 43   # contrabass


def peak_detection(wavenumbers, transmittances, min_depth: float = 20):
    """Detect IR absorption peaks as local minima in transmittance.

    IR absorption peaks appear as dips (local minima) in a transmittance
    spectrum. A peak is accepted if its absorption depth exceeds min_depth.

    Args:
        wavenumbers: list of wavenumbers in cm⁻¹
        transmittances: list of transmittance values in % (0–100)
        min_depth: minimum absorption depth (100 - transmittance) to be
                   counted as a peak. Default 20 filters out baseline noise.
    Returns:
        list of (wavenumber, depth) tuples, sorted by wavenumber descending
    """
    peaks = []
    for i in range(1, len(transmittances) - 1):
        is_local_min = (transmittances[i] < transmittances[i - 1] and
                        transmittances[i] < transmittances[i + 1])
        depth = 100 - transmittances[i]
        if is_local_min and depth > min_depth:
            peaks.append((wavenumbers[i], depth))

    peaks.sort(key=lambda p: p[0], reverse=True)
    return peaks


def wavenumber_to_midi(wavenumber: float) -> int:
    """Map wavenumber (cm⁻¹) to MIDI pitch using linear scaling.

    Higher wavenumber = faster vibration = higher pitch.
    Range: 400 cm⁻¹ → MIDI 40 (E2), 4000 cm⁻¹ → MIDI 84 (C6).

    Args:
        wavenumber: IR wavenumber in cm⁻¹
    Returns:
        MIDI pitch number (40–84)
    """
    ratio = (wavenumber - 400) / (4000 - 400)
    pitch = 40 + ratio * 44
    return int(round(max(40, min(84, pitch))))


def depth_to_velocity(depth: float) -> int:
    """Map absorption depth (0–100%) to MIDI velocity (45–110).

    Stronger absorption → louder note.

    Args:
        depth: absorption depth in % (= 100 - transmittance)
    Returns:
        MIDI velocity (45–110)
    """
    return int(45 + (depth / 100) * 65)


def depth_to_duration(depth: float) -> float:
    """Map absorption depth (0–100%) to note duration in beats (0.25–1.0).

    Stronger absorption = more prominent spectral feature = longer note,
    so the listener has time to register it.

    Args:
        depth: absorption depth in %
    Returns:
        note duration in beats
    """
    return round(0.25 + (depth / 100) * 0.75, 3)


def spectrum_to_midi(compound_or_cas, output_path: str, data=None, min_depth: float = 20):
    """Convert a molecule's IR spectrum into a MIDI file.

    Pipeline:
        compound → NIST IR spectrum → peak detection → wavenumber-to-pitch
        mapping → MIDI file. Tempo and instrument are set by molecular weight.

    Args:
        compound_or_cas: CAS string or already-fetched NistCompound object.
        output_path: file path to write the .mid file
        data: optional pre-computed (wavenumbers, transmittances) tuple.
              Pass this to avoid fetching the IR spectrum from NIST again.
        min_depth: minimum absorption depth (%) to include as a note (default 20)

    Example:
        >>> spectrum_to_midi('67-64-1', 'acetone.mid')
    """
    if isinstance(compound_or_cas, str):
        compound = nist.get_compound(compound_or_cas)
    else:
        compound = compound_or_cas

    if data is not None:
        wavenumbers, transmittances = data
    else:
        wavenumbers, transmittances = extract_spectrum_data(compound)

    bpm        = mw_to_bpm(compound.mol_weight)
    instrument = mw_to_instrument(compound.mol_weight)
    peaks      = peak_detection(wavenumbers, transmittances, min_depth=min_depth)

    midi = MIDIFile(1)
    midi.addTempo(0, 0, bpm)
    midi.addProgramChange(0, 0, 0, instrument)

    current_beat = 0.0
    for wavenumber, depth in peaks:
        pitch    = wavenumber_to_midi(wavenumber)
        velocity = depth_to_velocity(depth)
        duration = depth_to_duration(depth)
        midi.addNote(0, 0, pitch, current_beat, duration, velocity)
        current_beat += duration

    with open(output_path, "wb") as f:
        midi.writeFile(f)
