"""Tests for the ``utils`` module."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "chemecho"))

from utils import (
    _detect_input_type,
    _smiles_from,
    draw_molecule,
    molecule_3d_html,
    nist_compound_from,
    resolve_molecule,
)


def test_draw_molecule():
    # test the size of the image
    assert draw_molecule("CCO").size == (400, 300)
    # test with a wrong sequence, should return None
    assert draw_molecule("abc") == None
    # test with an empty sequence, should return None
    assert draw_molecule("") == None


class TestDetectInputType:

    # checks cas input
    def test_cas_number_detected(self):
        assert _detect_input_type("58-08-2") == "cas"

    # checks smiles input with equals (double bonds)
    def test_smiles_detected_by_equals(self):
        assert _detect_input_type("CC(=O)C") == "smiles"

    # checks smiles input with hash (triple bond)
    def test_smiles_detected_by_hash(self):
        assert _detect_input_type("CC#N") == "smiles"

    # check formula input
    def test_formula_detected(self):
        assert _detect_input_type("C6H12O6") == "formula"

    # check name input
    def test_name_detected(self):
        assert _detect_input_type("caffeine") == "name"


class TestSmilesFrom:

    # checks that when requesting in PubChem it prefers the isomeric SMILES which contains more info
    def test_isomeric_smiles_preferred(self):
        props = {"IsomericSMILES": "CC(=O)C", "CanonicalSMILES": "CCC"}
        assert _smiles_from(props) == "CC(=O)C"

    # checks that the canonical SMILES is used when isomeric not present
    def test_canonical_smiles_fallback(self):
        props = {"CanonicalSMILES": "CCC"}
        assert _smiles_from(props) == "CCC"

    # checks that if no SMILES exists it returns None
    def test_empty_dict_returns_none(self):
        assert _smiles_from({}) is None


class TestResolveMolecule:

    # checks that all info exists and is returned when entering name
    def test_name_query_works(self):
        result = resolve_molecule("caffeine")
        assert result["smiles"] is not None
        assert result["cas"] is not None
        assert result["name"] is not None

    # checks that all info exists and is returned when entering cas
    def test_cas_query_works(self):
        result = resolve_molecule("58-08-2")
        assert result["smiles"] is not None
        assert result["cas"] is not None
        assert result["name"] is not None

    # checks that all info exists and is returned when entering smiles
    def test_smiles_query_works(self):
        result = resolve_molecule("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        assert result["smiles"] is not None
        assert result["cas"] is not None
        assert result["name"] is not None

    # checks that all info exists and is returned when entering formula
    def test_formula_query_works(self):
        result = resolve_molecule(" C8H10N4O2")
        assert result["smiles"] is not None
        assert result["cas"] is not None
        assert result["name"] is not None

    # checks that error (PubChem can't find compound) raised on invalid entry
    def test_pubchem_failure_raises_value_error(self):
        with pytest.raises(ValueError, match="PubChem could not find"):
            resolve_molecule("notacompound12345xyz")


class TestNistCompoundFrom:

    # checks that NIST search works
    def test_cas_found_returns_compound(self):
        cas, compound = nist_compound_from("67-64-1", "acetone")
        assert cas == "67-64-1"
        assert compound is not None

    # checks that NIST search gives error when invalid search
    def test_no_nist_results_raises_value_error(self):
        with pytest.raises(ValueError, match="No IR data found"):
            nist_compound_from(None, "notacompound12345xyz")


class TestMolecule3dHtml:

    # checks that the html string returned exists (is not empty, len > 0) and is a string type
    def test_valid_smiles_returns_html_string(self):
        html = molecule_3d_html("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        assert isinstance(html, str)
        assert len(html) > 0

    # checks that None is returned when invalid SMILES
    def test_invalid_smiles_returns_none(self):
        assert molecule_3d_html("not_a_smiles!!!") is None
