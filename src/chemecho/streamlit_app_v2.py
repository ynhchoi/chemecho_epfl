import streamlit as st
import streamlit.components.v1 as components
import nistchempy as nist
import pandas as pd
import os
import io
import base64

from get_spectrum import extract_spectrum_data, ir_graph, good_spectrum_count
from musification import molecular_music
from utils import get_smiles, draw_molecule, name_to_cas




# ── UI ────────────────────────────────────────────────────────────────────────

st.title("ChemEcho 🎶")
st.write("Enter a molecule name to see its structure, IR spectrum, and hear it as music.")

name_or_cas = st.selectbox(
    "Research by:",
    options=["CAS (recommended)", "Name"]
)
if name_or_cas == "CAS (recommended)":
    cas = st.text_input("**Molecule CAS**", "58-08-2")

else :
    name = st.text_input("**Molecule name**", "acetone")
    cas = name_to_cas(name)

            

compound = nist.get_compound(cas)
#data = extract_spectrum_data(compound)
spectrum_count, good_index_list = good_spectrum_count(compound)

st.sidebar.markdown('# Options')
tempo = st.sidebar.slider("Spectrum music tempo", 20, 200, 120)
if spectrum_count == 1:
    st.caption(f"Only one spectrum available for {compound.name}")
    spectrum_index = good_index_list[0]
    data = extract_spectrum_data(compound, spectrum_index)
else:
    index = st.sidebar.slider("Number of spectrum that can be translated", 1, spectrum_count, 1)
    spectrum_index = good_index_list[index-1]
    data = extract_spectrum_data(compound, spectrum_index)


st.subheader(compound.name)
st.caption(f"CAS: {cas}")

col1, col2 = st.columns(2)

# Molecular structure
with col1:
    st.markdown("**Molecular structure**")
    smiles = get_smiles(cas)
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
    fig, ax = ir_graph(data, compound.name)
    st.pyplot(fig)
    buffer = io.StringIO()
    fig.savefig(buffer, format="svg")
    svg = buffer.getvalue()

    st.download_button(
        label="Download spectrum as SVG",
        data=svg,
        file_name=f"{compound.name}_spectrum.csv",
        mime="image/svg+xml")


# MIDI generation
try:
    midi_path = molecular_music(data, compound, tempo)

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
""", height=150)

st.download_button(
    label="Download MIDI",
    data=midi_bytes,
    file_name=f"{compound.name}.mid",
    mime="audio/midi")
