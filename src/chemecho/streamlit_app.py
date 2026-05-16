import streamlit as st
import streamlit.components.v1 as components
import nistchempy as nist
import pandas as pd
import requests
import tempfile
import os
import base64
from rdkit import Chem
from rdkit.Chem import Draw

from get_spectrum import extract_spectrum_data, ir_graph, from_df_to_csv
import importlib.util, pathlib
_mus_spec = importlib.util.spec_from_file_location(
    "musification",
    pathlib.Path(__file__).parent / "musification_2.0.py"
)
_mus_mod = importlib.util.module_from_spec(_mus_spec)
_mus_spec.loader.exec_module(_mus_mod)
molecular_music = _mus_mod.molecular_music


def name_to_cas(name: str) -> str:
    """Search NIST by compound name and return CAS number.

    NIST compound IDs follow the format C{CAS_without_dashes}, e.g.
    acetone (67-64-1) has ID 'C67641'. We reconstruct the CAS as:
        numeric[:-3] - numeric[-3:-1] - numeric[-1]
    """
    results = nist.run_search(name, 'name', cIR=True)
    if not results.compounds:
        raise ValueError(f"No compound with IR data found for '{name}' on NIST.")
    nist_id = results.compounds[0].ID       # e.g. 'C67641&Units=SI'
    numeric = nist_id.lstrip('C').split('&')[0]   # -> '67641'
    cas = f"{numeric[:-3]}-{numeric[-3:-1]}-{numeric[-1]}"
    return cas


def get_smiles(name: str) -> str | None:
    """Fetch canonical SMILES from PubChem by compound name."""
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{name}/property/IsomericSMILES/JSON"
    )
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            props = resp.json()['PropertyTable']['Properties'][0]
            return props.get('IsomericSMILES') or props.get('SMILES')
    except Exception:
        pass
    return None


def draw_molecule(smiles: str):
    """Return a PIL image of the molecule from its SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Draw.MolToImage(mol, size=(400, 300))


# ── UI ────────────────────────────────────────────────────────────────────────

st.title("Chem Echo 🎶")
st.write("Enter a molecule name to see its structure, IR spectrum, and hear it as music.")

name = st.text_input("Molecule name", placeholder="e.g. acetone, ethanol, caffeine")

if st.button("Generate") and name:
    with st.spinner("Fetching data..."):
        try:
            # Name → CAS
            cas = name_to_cas(name)
            compound = nist.get_compound(cas)
            st.subheader(compound.name)
            st.caption(f"CAS: {cas}")

            col1, col2 = st.columns(2)

            # Molecular structure
            with col1:
                st.markdown("**Molecular structure**")
                smiles = get_smiles(name)
                if smiles:
                    img = draw_molecule(smiles)
                    if img:
                        st.image(img)
                    st.caption(f"SMILES: `{smiles}`")
                else:
                    st.warning("Structure not available from PubChem.")

            # IR spectrum
            with col2:
                st.markdown("**IR spectrum**")
                data = extract_spectrum_data(compound)
                fig, ax = ir_graph(data, compound.name)
                st.pyplot(fig)

            # CSV download
            df = pd.DataFrame({"Wavenumber": data[0], "Transmittance": data[1]})
            st.download_button(
                label="Download spectrum as CSV",
                data=from_df_to_csv(df),
                file_name=f"{cas}_spectrum.csv",
                mime="text/csv"
            )

            # MIDI generation
            with tempfile.NamedTemporaryFile(suffix=".mid", delete=False) as tmp:
                midi_path = tmp.name

            try:
                molecular_music(compound, midi_path, data=data)
                with open(midi_path, "rb") as f:
                    midi_bytes = f.read()
            finally:
                if os.path.exists(midi_path):
                    os.remove(midi_path)

            st.markdown("**Listen**")
            midi_b64 = base64.b64encode(midi_bytes).decode()
            components.html(f"""
                <script src="https://cdn.jsdelivr.net/npm/tone@14/build/Tone.js"></script>
                <script src="https://cdn.jsdelivr.net/npm/@magenta/music@1.23.1/es6/core.js"></script>
                <script src="https://cdn.jsdelivr.net/npm/html-midi-player@1.5.0"></script>
                <midi-player
                    src="data:audio/midi;base64,{midi_b64}"
                    sound-font
                    style="width:100%;margin-top:4px">
                </midi-player>
            """, height=120)

            st.download_button(
                label="Download MIDI",
                data=midi_bytes,
                file_name=f"{compound.name}.mid",
                mime="audio/midi"
            )

        except Exception as e:
            st.error(f"Something went wrong: {e}")
