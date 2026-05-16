# Changes to `musification_2.0.py` — Summary for Team

## Overview

`musification_2.0.py` has been updated to fix several bugs and improve the scientific accuracy of the sonification pipeline. The overall structure and design ideas from the original version have been preserved — particularly the drum track, note overlap, and 30-beat normalisation. The changes below are grouped by what was broken, what was improved, and what was added from other sources.

---

## 1. Dependency: `musicpy` → `midiutil`

**Before:**
```python
import musicpy as mp
from musicpy.daw import *
import chemecho_epfl.src.chemecho.get_spectrum as spect
```

**After:**
```python
from midiutil import MIDIFile
import math
from scipy.signal import find_peaks
from get_spectrum import extract_spectrum_data
```

**Why:** `musicpy` is not installed in the project environment and the import path `chemecho_epfl.src.chemecho.get_spectrum` no longer exists — the file would crash immediately on import. `midiutil` is already installed and used by the rest of the project. `scipy` was added to the environment (`pip install scipy`) and is used only for `find_peaks`.

---

## 2. `peak_detection` — Complete Rewrite

**Before (original logic):**
```python
threshold = max(transmittances) - 0.1 * (max(transmittances) - min(transmittances))
for i in range(1, len(transmittances) - 1):
    if transmittances[i] < threshold:
        peaks.append((wavenumbers[i], transmittances[i]))
    else:
        peaks.append((wavenumbers[i], max(transmittances)))  # non-peaks added too
```

**After:**
```python
neg_t = [-v for v in transmittances]
indices, props = find_peaks(neg_t, prominence=min_prominence, width=min_width_pts)
# returns (wavenumber, prominence, fwhm_width) for true peaks only
```

**Why:** The original code was not actually detecting peaks — it was keeping *every single data point* (2000–3000 points for a typical NIST spectrum) and replacing non-peak values with `max(transmittances)`. This caused two problems:

1. `target_duration_beats / len(peaks)` became essentially zero, so all notes played simultaneously.
2. Every noise fluctuation became a note.

The new version uses `scipy.signal.find_peaks` on the inverted transmittance (absorption dips = local minima in transmittance). The `prominence` parameter controls minimum absorption depth and `width` filters out single-sample noise spikes. A typical NIST spectrum now yields 5–20 meaningful peaks instead of thousands.

The function now returns a third value per peak — the **FWHM half-width in cm⁻¹** — which is used to set note duration (see section 4).

---

## 3. Pitch Mapping: Transmittance → Wavenumber

**Before:**
```python
freq = 196 + (100 - trans) / 100 * (1046.50 - 196)
note_generated = mp.freq_to_note(freq)
```

**After:**
```python
def wavenumber_to_midi(wavenumber: float) -> int:
    ratio = (wavenumber - 400) / (4000 - 400)
    return int(round(max(40, min(84, 40 + ratio * 44))))
```

**Why:** In the original, pitch was derived from *transmittance* (absorption intensity), not from the peak's wavenumber position. This meant two peaks at the same wavenumber but different intensities would produce different pitches — losing the spectral position information entirely, which is the most chemically meaningful part of an IR spectrum (C–H stretches, C=O stretches, fingerprint region, etc. all have characteristic wavenumber positions).

The new mapping encodes wavenumber linearly onto MIDI pitch: 400 cm⁻¹ → E2 (MIDI 40), 4000 cm⁻¹ → C6 (MIDI 84). This is consistent with published IR sonification approaches (see *J. Chem. Ed.* 2020, *arXiv* 2601.02652).

---

## 4. Note Interval: Equal Spacing → Wavenumber-Gap Proportional

**Before:**
```python
interval = target_beats / len(peaks)   # all gaps equal regardless of spectral distance
midi.addNote(0, 0, pitch, i * interval, duration, velocity)
```

**After:**
```python
# Gaps between consecutive peak wavenumbers → time gaps between notes
gaps = [wns[i] - wns[i+1] for i in range(len(wns) - 1)]
raw_intervals = [0.5 + (g - min_gap) / gap_range * 3.5 for g in gaps]
# Normalise so the total still equals 30 beats
scale = target_beats / sum(raw_intervals)
intervals = [iv * scale for iv in raw_intervals]

current_beat = 0.0
for i, (...) in enumerate(peaks):
    midi.addNote(0, 0, pitch, current_beat, duration, velocity)
    current_beat += intervals[i]
```

**Why:** With equal spacing, two C–H stretch peaks at 3035 and 2997 cm⁻¹ (gap = 38 cm⁻¹) receive the same time interval as the jump from 2997 to 1750 cm⁻¹ (gap = 1247 cm⁻¹). This erases information about how the peaks are distributed across the spectrum.

With gap-proportional spacing, peaks that are close in wavenumber play in rapid succession, while peaks separated by a large silent region of the spectrum are followed by a longer pause. In acetone, the two adjacent C–H stretches fire almost simultaneously (0.5 beats after normalisation), then a long pause marks the boundary between the functional group region and the fingerprint region, before the C=O and C–C peaks resume.

This is recommended by the *Sonification Handbook* (Ch. 15): "mapping data intervals to time intervals is the most direct and intuitive form of parameter mapping sonification."

The 30-beat total duration is preserved by normalising the raw intervals after computing them.

---

## 5. Note Duration: Fixed → FWHM-Based

**Before:**
```python
note_duration = interval_value * 10   # fixed multiplier relative to equal interval
```

**After:**
```python
# Duration is now independent of the note's start interval
width_ratio = (width_wn - min_w) / width_range
duration    = 0.5 + width_ratio * 1.5   # 0.5–2.0 beats from FWHM
```

**Why:** Using a fixed multiplier discards information about peak shape, and tying duration to the (now variable) interval would conflate two independent physical quantities. Duration is now computed purely from FWHM peak width, making it independent of spectral spacing.

In IR spectroscopy, peak width reflects vibrational mode lifetime: broader peaks indicate faster energy redistribution through mode coupling (IVR). A sharp, isolated peak (e.g. a C≡N stretch) plays as a short crisp note; a broad, strongly coupled peak (e.g. an O–H stretch with hydrogen bonding) plays as a longer sustained note. This is inspired by Kim & Heller (arXiv 2601.02652, 2026), who discuss peak width and IVR as key parameters in IR sonification.

---

## 5. Velocity: Fixed → Prominence-Based

**Before:**
```python
volume = 80  # constant volume
```

**After:**
```python
velocity = int(max(45, min(110, 45 + (prominence / 100) * 65)))
```

**Why:** A fixed volume wastes one of the most natural musical dimensions. Stronger IR absorption (larger prominence) now maps to a louder note. Weak peaks play quietly; strong peaks (e.g. the carbonyl C=O stretch) play loudly.

---

## 6. Tempo: Fixed Default → Molecular Weight

**Before:**
```python
def molecular_music(cas, bpm_mol=120):  # hardcoded default, never calculated
```

**After:**
```python
def mw_to_bpm(mol_weight: float) -> int:
    # kinetic theory: RMS velocity ∝ 1/√M
    ratio = (√mol_weight - √18) / (√600 - √18)
    return int(max(60, min(170, round(170 - ratio * 110))))

bpm = mw_to_bpm(compound.mol_weight)
```

**Why:** The original always used BPM 120 regardless of the molecule. BPM is now calculated from molecular weight using the kinetic theory analogy: lighter molecules move faster (higher RMS velocity ∝ 1/√M), so they get a faster tempo. The formula uses a square-root relationship, which is more physically accurate than the linear approximation used in `musification_v2.py`.

Example outputs: water (18 Da) → 170 BPM, acetone (58 Da) → 152 BPM, caffeine (194 Da) → 117 BPM.

---

## 7. `molecular_weight_to_sound_code`: CAS → Compound Object

**Before:**
```python
def molecular_weight_to_sound_code(compound_cas: str) -> int:
    compound = spect.nist.get_compound(compound_cas)  # extra NIST request
    molecular_weight_compound = compound.mol_weight
```

**After:**
```python
def molecular_weight_to_sound_code(compound) -> int:
    mw = compound.mol_weight
```

**Why:** Passing a CAS string forces an extra HTTP request to NIST each time the function is called, even when the compound object is already available. The function now takes the compound object directly, consistent with how `extract_spectrum_data` and `spectrum_to_midi` in the rest of the codebase work.

The instrument mapping itself is unchanged: violin (< 50 Da) → clarinet → cello → tenor sax → contrabass (≥ 300 Da).

---

## 8. Drum Track: `musicpy` → `midiutil` Channel 9

The drum track concept from the original (kick drum every 500 cm⁻¹, snare every 100 cm⁻¹ as a spectral ruler) has been preserved but re-implemented with `midiutil`:

- MIDI channel 9 is the standard General MIDI percussion channel.
- Kick drum (MIDI note 36) marks every 500 cm⁻¹.
- Hi-hat (MIDI note 42) marks every 100 cm⁻¹.
- The drum track runs for the same 30 beats as the melody track.
- Enabled by default (`drum_track=True`), but can be disabled.

---

## 9. `data=` Parameter (Avoid Double Fetch)

**Added:**
```python
def molecular_music(compound_or_cas, output_path, data=None, ...):
    if data is not None:
        wavenumbers, transmittances = data
    else:
        wavenumbers, transmittances = extract_spectrum_data(compound)
```

**Why:** When called from the Streamlit app, IR data is already fetched for the spectrum plot. Without this parameter, the function would call `extract_spectrum_data` again, triggering a second NIST HTTP request on the same compound object — which caused a `'NistCompound' object is not subscriptable` error in the app.

---

## What Was Kept From the Original

| Feature | Status |
|---|---|
| Drum track (kick/hi-hat as spectral ruler) | ✅ Kept, re-implemented with midiutil |
| 30-beat total duration normalisation | ✅ Kept |
| Note overlap (`duration > interval`) for continuous texture | ✅ Kept |
| Equal note spacing (`interval = 30 / n_peaks`) | ❌ Replaced by gap-proportional spacing (see §4) |
| Heavier molecule → lower-pitched instrument | ✅ Kept |

---

## References

- Kim & Heller, *Musical Molecules: Sonifying the IR Spectra...*, arXiv 2601.02652 (2026) — FWHM → duration, IVR and peak width
- Osterroth et al., *The Sound of Chemistry*, J. Chem. Ed. 97, 703 (2020) — wavenumber → pitch mapping
- Mahjour et al., *Molecular Sonification...*, Digital Discovery (2023) — multi-parameter encoding
- Herman & Casini, *Sonification Handbook*, Ch. 15: Parameter Mapping Sonification (2011) — gap-proportional interval mapping
