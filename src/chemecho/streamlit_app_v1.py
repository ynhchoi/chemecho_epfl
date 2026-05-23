import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import os
import io
import base64


from get_spectrum import extract_spectrum_data, ir_graph
from musification import molecular_music
from musification_fg import molecular_music_fg
from utils import resolve_molecule, nist_compound_from, draw_molecule, molecule_3d_html


# ── Sidebar ───────────────────────────────────────────────────────────────────

st.sidebar.markdown("# Options")
musification_version = st.sidebar.radio(
    "Musification",
    ["v1 (basic)", "v2 (functional groups)"],
    index=1,
)
if musification_version == "v1 (basic)":
    tempo = st.sidebar.slider("Tempo (BPM)", 20, 200, 120)


# ── UI ────────────────────────────────────────────────────────────────────────

st.title("Chem Echo 🎶")
st.write("Enter a molecule name, CAS number, molecular formula, or SMILES to see its structure, IR spectrum, and hear it as music.")

query = st.text_input("Molecule", placeholder="e.g. caffeine  |  67-64-1  |  C6H6  |  CC(=O)C")

if st.button("Generate") and query:
    with st.spinner("Fetching data..."):
        try:
            info = resolve_molecule(query)
            smiles = info['smiles']
            cas, compound = nist_compound_from(info['cas'], info['name'])
            extracted_data = extract_spectrum_data(compound)

            st.subheader(compound.name)
            st.caption(f"CAS: {cas}")

            col1, col2 = st.columns(2)

            # Molecular structure
            with col1:
                st.markdown("**Molecular structure**")
                if smiles:
                    tab_2d, tab_3d = st.tabs(["2D", "3D"])
                    with tab_2d:
                        img = draw_molecule(smiles)
                        if img:
                            st.image(img)
                    with tab_3d:
                        html_3d = molecule_3d_html(smiles)
                        if html_3d:
                            components.html(html_3d, height=400, scrolling=False)
                        else:
                            st.warning("3D structure could not be generated.")
                    st.caption(f"SMILES: `{smiles}`")
                else:
                    st.warning("Structure not available from PubChem.")

            # IR spectrum
            with col2:
                st.markdown("**IR spectrum**")
                fig, ax = ir_graph(extracted_data, compound.name)
                st.pyplot(fig)
                buffer = io.StringIO()
                fig.savefig(buffer, format="svg")
                svg = buffer.getvalue()

                st.download_button(
                    label="Download spectrum as SVG",
                    data=svg,
                    file_name=f"{compound.name}_spectrum.svg",
                    mime="image/svg+xml")

            # MIDI generation
            midi_path = None
            legend = None
            try:
                if musification_version == "v1 (basic)":
                    midi_path = molecular_music(extracted_data, compound, tempo)
                else:
                    if smiles:
                        midi_path, legend = molecular_music_fg(extracted_data, compound, smiles)
                    else:
                        st.warning("SMILES not available — falling back to v1.")
                        midi_path = molecular_music(extracted_data, compound, 120)

                with open(midi_path, "rb") as f:
                    midi_bytes = f.read()
            finally:
                if midi_path and os.path.exists(midi_path):
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
                mime="audio/midi"
            )

            # FG legend (v2 only)
            if legend:
                st.caption(f"BPM: {legend['bpm']}  ·  Carbon count: {legend['carbon_count']}")
                st.markdown("**Detected functional groups**")
                rows = [
                    {
                        "Functional group": fg,
                        "Region (cm⁻¹)": f"{info['region'][0]}–{info['region'][1]}",
                        "Peaks detected": info['n_peaks'],
                    }
                    for fg, info in legend['fgs'].items()
                ]
                st.table(pd.DataFrame(rows))

        except Exception as e:
            st.error(f"Something went wrong: {e}")
