import os
import pytest
import sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import nistchempy as nist
import matplotlib.pyplot as plt


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "chemecho"))


from get_spectrum import extract_spectrum_data, ir_graph
from utils import (
    _detect_input_type,
    _smiles_from,
    draw_molecule,
    molecule_3d_html,
    nist_compound_from,
    resolve_molecule,
)
from musification_fg import (  # noqa: E402
    FG_CATALOG,
    _scale_accent_volume,
    assign_peaks_to_fgs,
    count_carbons,
    detect_functional_groups,
    find_absorption_peaks,
    molecular_music_fg,
    molecular_weight_to_bpm,
    noise_reduction
)


# TESTS FOR FUNCTIONS FROM get_spectrum FILE

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
    assert pytest.approx(ax.get_xlim(),abs=10) == (3966, 450)

    plt.close(fig)


# TESTS FOR FUNCTIONS FROM utils FILE

def test_draw_molecule():
    #test the size of the image
    assert draw_molecule("CCO").size == (400, 300)
    #test with a wrong sequence, should return "None"
    assert draw_molecule("abc") == None
    #test with an empty sequence, should return "None"
    assert draw_molecule("") == None

class TestDetectInputType:

    #checks cas input
    def test_cas_number_detected(self):
        assert _detect_input_type("58-08-2") == "cas"

    #checks smiles input with equals (double bonds)
    def test_smiles_detected_by_equals(self):
        assert _detect_input_type("CC(=O)C") == "smiles"

    #checks smules input with hash (triple bond)
    def test_smiles_detected_by_hash(self):
        assert _detect_input_type("CC#N") == "smiles"

    #check formula input
    def test_formula_detected(self):
        assert _detect_input_type("C6H12O6") == "formula"

    #check name input
    def test_name_detected(self):
        assert _detect_input_type("caffeine") == "name"
   
class TestSmilesFrom:

    #checks that when requesting in Pubchem it prefers the isomeric smiles which contains more info
    def test_isomeric_smiles_preferred(self):
        props = {"IsomericSMILES": "CC(=O)C", "CanonicalSMILES": "CCC"}
        assert _smiles_from(props) == "CC(=O)C"

    #checks that the canonical smiles is used when ismoeric not present
    def test_canonical_smiles_fallback(self):
        props = {"CanonicalSMILES": "CCC"}
        assert _smiles_from(props) == "CCC"

    #checks that if no smiles exists it returns None
    def test_empty_dict_returns_none(self):
        assert _smiles_from({}) is None

class TestResolveMolecule:

    #checks that all info exists and is returned when entering name
    def test_name_query_works(self):
        result = resolve_molecule("caffeine")
        assert result["smiles"] is not None
        assert result["cas"] is not None
        assert result["name"] is not None

    #checks that all info exists and is returned when entering cas
    def test_name_query_works(self):
        result = resolve_molecule("58-08-2")
        assert result["smiles"] is not None
        assert result["cas"] is not None
        assert result["name"] is not None

    #checks that all info exists and is returned when entering smiles
    def test_name_query_works(self):
        result = resolve_molecule("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        assert result["smiles"] is not None
        assert result["cas"] is not None
        assert result["name"] is not None
    
    #checks that all info exists and is returned when entering fromula
    def test_name_query_works(self):
        result = resolve_molecule(" C8H10N4O2")
        assert result["smiles"] is not None
        assert result["cas"] is not None
        assert result["name"] is not None
    
    #checks that error (Pubchem cant find compound) if invalid entry
    def test_pubchem_failure_raises_value_error(self):
        with pytest.raises(ValueError, match="PubChem could not find"):
            resolve_molecule("notacompound12345xyz")

class TestNistCompoundFrom:

    #checks that NIST search works
    def test_cas_found_returns_compound(self):
        cas, compound = nist_compound_from("67-64-1", "acetone")
        assert cas == "67-64-1"
        assert compound is not None

    #checks that NSIT search gives error when invalid search
    def test_no_nist_results_raises_value_error(self):
        with pytest.raises(ValueError, match="No IR data found"):
            nist_compound_from(None, "notacompound12345xyz")

class TestMolecule3dHtml:

    #checks that html string returned exists (is not nothing , len>0) and is a string type
    def test_valid_smiles_returns_html_string(self):
        html = molecule_3d_html("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        assert isinstance(html, str)
        assert len(html) > 0

    #checks that None is returned when invalid smiles 
    def test_invalid_smiles_returns_none(self):
        assert molecule_3d_html("not_a_smiles!!!") is None

# TESTS FOR FUNCTIONS FROM musification_fg FILEQ


# ─── synthetic spectrum helper ────────────────────────────────────────────────

def make_synthetic_spectrum(peak_cm1, depths, width=20, wn_range=(500, 4000), n=1000):
    """Build a (wavenumbers, transmittances) pair with Gaussian dips at given cm-1.

    Wavenumbers run from high to low to match NIST convention.
    """
    wavenumbers = np.linspace(wn_range[1], wn_range[0], n)
    transmittance = np.full(n, 100.0)
    for pos, depth in zip(peak_cm1, depths):
        transmittance -= depth * np.exp(-((wavenumbers - pos) ** 2) / (2 * width ** 2))
    return wavenumbers.tolist(), transmittance.tolist()


# ─── detect_functional_groups ─────────────────────────────────────────────────

class TestDetectFunctionalGroups:
    def test_empty_string_returns_empty(self):
        assert detect_functional_groups("") == []

    def test_invalid_smiles_returns_empty(self):
        assert detect_functional_groups("not_a_smiles_!!!") == []

    def test_acetone_has_carbonyl_and_csp3(self):
        fgs = detect_functional_groups("CC(=O)C")
        assert "C=O (carbonyl)" in fgs
        assert "C-H (sp3)" in fgs

    def test_benzene_has_aromatic_ring_and_aromatic_ch(self):
        fgs = detect_functional_groups("c1ccccc1")
        assert "aromatic ring" in fgs
        assert "C-H (aromatic)" in fgs

    def test_methanol_has_alcohol_and_co(self):
        fgs = detect_functional_groups("CO")
        assert "O-H (alcohol)" in fgs
        assert "C-O" in fgs

    def test_acetonitrile_has_nitrile(self):
        assert "C#N (nitrile)" in detect_functional_groups("CC#N")

    def test_acetic_acid_has_carboxylic_oh(self):
        fgs = detect_functional_groups("CC(=O)O")
        assert "O-H (carboxylic)" in fgs
        assert "C=O (carbonyl)" in fgs

    def test_isolated_water_atom_matches_nothing(self):
        # The catalog's [OX2H] requires the O to have 2 graph connections,
        # so a lone water O (no heavy neighbors) is correctly NOT matched as alcohol.
        assert detect_functional_groups("O") == []

    def test_alkyne_detected(self):
        assert "C#C (alkyne)" in detect_functional_groups("CC#CC")

    def test_nitrobenzene_has_nitro(self):
        # RDKit canonicalizes nitro as the zwitterionic [N+](=O)[O-] form.
        assert "N=O (nitro)" in detect_functional_groups("c1ccccc1[N+](=O)[O-]")


# ─── count carbons ────────────────────────────────────────────────────────────

class TestCountCarbons:
    @pytest.mark.parametrize("smiles,expected", [
        ("", 0),
        ("invalid_xyz!!", 0),
        ("C", 1),                    # methane
        ("CC", 2),                   # ethane
        ("CC(=O)C", 3),              # acetone
        ("c1ccccc1", 6),             # benzene
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 8),  # caffeine
        ("O", 0),                    # water (no C)
        ("N#N", 0),                  # N2 gas
    ])
    def test_known_carbon_counts(self, smiles, expected):
        assert count_carbons(smiles) == expected


# ─── find absorption peaks ────────────────────────────────────────────────────

class TestFindAbsorptionPeaks:
    def test_flat_spectrum_returns_empty(self):
        wn = list(np.linspace(4000, 500, 500))
        trans = [100.0] * 500
        assert find_absorption_peaks(wn, trans) == []

    def test_single_peak_detected(self):
        wn, trans = make_synthetic_spectrum(peak_cm1=[1700], depths=[60])
        peaks = find_absorption_peaks(wn, trans)
        assert len(peaks) == 1
        assert abs(peaks[0]["wavenumber"] - 1700) < 30

    def test_multiple_peaks_detected(self):
        wn, trans = make_synthetic_spectrum(
            peak_cm1=[3300, 2950, 1700, 1100], depths=[50, 40, 70, 30]
        )
        peaks = find_absorption_peaks(wn, trans)
        wn_found = sorted(p["wavenumber"] for p in peaks)
        for target in [1100, 1700, 2950, 3300]:
            assert any(abs(w - target) < 30 for w in wn_found), f"peak near {target} not found"

    def test_low_prominence_peak_filtered_out(self):
        # Prominence threshold is relative to the spectrum's range. With one big
        # peak (depth 50) and one tiny peak (depth 1, which is 2% of range), the
        # tiny one should be filtered out by the default 5% threshold.
        wn, trans = make_synthetic_spectrum(
            peak_cm1=[1700, 2900], depths=[50, 1]
        )
        peaks = find_absorption_peaks(wn, trans)
        wavenumbers_found = [p["wavenumber"] for p in peaks]
        assert any(abs(w - 1700) < 30 for w in wavenumbers_found), "big peak should be found"
        assert not any(abs(w - 2900) < 30 for w in wavenumbers_found), "tiny peak should be filtered"

    def test_peaks_have_required_keys(self):
        wn, trans = make_synthetic_spectrum(peak_cm1=[1700], depths=[60])
        peaks = find_absorption_peaks(wn, trans)
        assert peaks, "expected at least one peak"
        for key in ("slot", "wavenumber", "transmittance", "prominence"):
            assert key in peaks[0]

    def test_slot_indexes_music_timeline(self):
        # slot = original_index - 1 (because noise_reduction skips index 0 and N-1)
        wn, trans = make_synthetic_spectrum(peak_cm1=[1700], depths=[60])
        peaks = find_absorption_peaks(wn, trans)
        for p in peaks:
            assert p["slot"] == wn.index(p["wavenumber"]) - 1


# ─── assign peaks to functional groups ──────────────────────────────────────────────────────

class TestAssignPeaksToFgs:
    @staticmethod
    def _peak(wn, slot=0, trans=50.0, prom=30.0):
        return {
            "wavenumber": wn,
            "slot": slot,
            "transmittance": trans,
            "prominence": prom,
        }

    def test_empty_peaks_returns_empty_dict(self):
        assert assign_peaks_to_fgs([], ["C=O (carbonyl)"]) == {}

    def test_empty_fgs_returns_empty_dict(self):
        assert assign_peaks_to_fgs([self._peak(1700)], []) == {}

    def test_peak_in_single_region_assigned(self):
        result = assign_peaks_to_fgs([self._peak(1700)], ["C=O (carbonyl)"])
        assert "C=O (carbonyl)" in result
        assert len(result["C=O (carbonyl)"]) == 1

    def test_peak_outside_all_regions_dropped(self):
        # 500 cm-1 is below every FG region — should be dropped
        result = assign_peaks_to_fgs([self._peak(500)], list(FG_CATALOG.keys()))
        assert result == {}

    def test_narrowest_region_wins_disambiguation(self):
        # 2900 cm-1 falls inside O-H (carboxylic) [2500–3300] AND C-H (sp3) [2840–2985].
        # C-H sp3 is narrower (145 cm-1) than O-H carb (800 cm-1), so it should win.
        peaks = [self._peak(2900)]
        present = ["O-H (carboxylic)", "C-H (sp3)"]
        result = assign_peaks_to_fgs(peaks, present)
        assert list(result.keys()) == ["C-H (sp3)"]

    def test_multiple_peaks_assigned_to_same_fg(self):
        peaks = [self._peak(1700), self._peak(1750)]
        result = assign_peaks_to_fgs(peaks, ["C=O (carbonyl)"])
        assert len(result["C=O (carbonyl)"]) == 2


# ─── scale accent volume ─────────────────────────────────────────────────────

class TestScaleAccentVolume:
    def test_max_prom_zero_returns_hi(self):
        assert _scale_accent_volume(0.0, 0.0) == 127

    def test_zero_prominence_returns_lo(self):
        assert _scale_accent_volume(0.0, 10.0, lo=60, hi=127) == 60

    def test_full_prominence_returns_hi(self):
        assert _scale_accent_volume(10.0, 10.0, lo=60, hi=127) == 127

    def test_half_prominence_returns_midpoint(self):
        v = _scale_accent_volume(5.0, 10.0, lo=60, hi=127)
        assert v == int(60 + (127 - 60) * 0.5)

    def test_custom_range(self):
        assert _scale_accent_volume(1.0, 1.0, lo=0, hi=100) == 100
        assert _scale_accent_volume(0.0, 1.0, lo=0, hi=100) == 0


# ─── molecular weight to bpm ──────────────────────────────────────────────────

class TestMolecularWeightToBpm:
    def test_light_molecule_returns_bpm_max(self):
        # M_fast reference is 16, so MW=16 should map to upper bound
        compound = SimpleNamespace(mol_weight=16.0)
        assert molecular_weight_to_bpm(compound, bpm_min=60, bpm_max=120) == 120

    def test_heavy_molecule_returns_bpm_min(self):
        compound = SimpleNamespace(mol_weight=500.0)
        assert molecular_weight_to_bpm(compound, bpm_min=60, bpm_max=120) == 60

    def test_extra_heavy_clamped_to_bpm_min(self):
        compound = SimpleNamespace(mol_weight=10_000.0)
        assert molecular_weight_to_bpm(compound, bpm_min=60, bpm_max=120) == 60

    def test_extra_light_clamped_to_bpm_max(self):
        compound = SimpleNamespace(mol_weight=1.0)
        assert molecular_weight_to_bpm(compound, bpm_min=60, bpm_max=120) == 120

    def test_intermediate_mass_is_between_bounds(self):
        compound = SimpleNamespace(mol_weight=100.0)
        bpm = molecular_weight_to_bpm(compound, bpm_min=60, bpm_max=120)
        assert 60 < bpm < 120

    def test_custom_bounds(self):
        compound = SimpleNamespace(mol_weight=16.0)
        assert molecular_weight_to_bpm(compound, bpm_min=80, bpm_max=160) == 160


# ─── molecular music fg (integration with mocked compound) ────────────────────

class TestMolecularMusicFg:
    def test_writes_midi_and_returns_legend(self, tmp_path, monkeypatch):
        # Force MIDI to be written into tmp_path so the test cleans up automatically.
        monkeypatch.chdir(tmp_path)

        compound = SimpleNamespace(name="testmol", mol_weight=100.0)
        smiles = "CC(=O)C"  # acetone — has C=O and C-H sp3
        # Synthetic spectrum with peaks at carbonyl (~1715) and CH-sp3 (~2900) regions.
        wn, trans = make_synthetic_spectrum(peak_cm1=[1715, 2900], depths=[70, 50])

        filename, legend = molecular_music_fg((wn, trans), compound, smiles)

        assert (tmp_path / filename).exists()
        assert legend["carbon_count"] == 3
        assert 60 <= legend["bpm"] <= 120
        assert "C=O (carbonyl)" in legend["fgs"]

    def test_invalid_smiles_yields_no_fgs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        compound = SimpleNamespace(name="x", mol_weight=100.0)
        wn, trans = make_synthetic_spectrum(peak_cm1=[1700], depths=[60])
        _, legend = molecular_music_fg((wn, trans), compound, "not_smiles")
        assert legend["fgs"] == {}
        assert legend["carbon_count"] == 0

# tests for noise reduction

class TestNoiseReduction:

    #checks that the output of noise_reduction is a list of tuples
    def test_returns_list_of_tuples(self):
        wn, trans = make_synthetic_spectrum(peak_cm1=[1700], depths=[60])
        result = noise_reduction(wn, trans)
        assert isinstance(result, list)
        assert all(isinstance(item, tuple) and len(item) == 2 for item in result)

    # checks that the noise reduction output doesnt take the first and last data points which can sometimes cause problems
    def test_skips_first_and_last_sample(self):
        wn, trans = make_synthetic_spectrum(peak_cm1=[1700], depths=[60])
        result = noise_reduction(wn, trans)
        assert len(result) == len(wn) - 2

    #checks that the transmittances to be ignored (the very high ones, essentially noise) are actually deleted and set to max transmittance
    def test_noise_replaced_by_max_transmittance(self):
        wn, trans = make_synthetic_spectrum(peak_cm1=[1700], depths=[60])
        max_trans = max(trans)
        threshold = max_trans - 0.1 * (max_trans - min(trans))
        result = noise_reduction(wn, trans)
        for i, (w, t) in enumerate(result):
            original_trans = trans[i + 1]  # +1 because first sample is skipped
            if original_trans >= threshold:
                assert t == max_trans

    # checks that the peaks that should still remain actually stay
    def test_real_peaks_preserved(self):
        wn, trans = make_synthetic_spectrum(peak_cm1=[1700], depths=[60])
        max_trans = max(trans)
        threshold = max_trans - 0.1 * (max_trans - min(trans))
        result = noise_reduction(wn, trans)
        for i, (w, t) in enumerate(result):
            original_trans = trans[i + 1]
            if original_trans < threshold:
                assert t == original_trans

    # checks that the wavenumbers data is unchanged
    def test_wavenumbers_preserved(self):
        wn, trans = make_synthetic_spectrum(peak_cm1=[1700], depths=[60])
        result = noise_reduction(wn, trans)
        result_wn = [w for w, t in result]
        assert result_wn == wn[1:-1]

    # checks the case of a flat spectrum with no actual peaks so all transmittances are the max (which is the same as all tranmittances)
    def test_flat_spectrum_all_replaced_by_max(self):
        wn = list(np.linspace(4000, 500, 100))
        trans = [80.0] * 100
        result = noise_reduction(wn, trans)
        assert all(t == 80.0 for _, t in result)


    # checks the case of a spectrum with only 1 peak 
    def test_single_deep_peak_survives(self):
        wn, trans = make_synthetic_spectrum(peak_cm1=[1700], depths=[80], width=5)
        max_trans = max(trans)
        result = noise_reduction(wn, trans)
        low_trans_points = [(w, t) for w, t in result if t < max_trans]
        assert len(low_trans_points) > 0
        assert all(abs(w - 1700) < 50 for w, t in low_trans_points)
