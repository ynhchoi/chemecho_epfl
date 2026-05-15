import streamlit as st
import nistchempy as nist
import tempfile
import os

from chemecho_epfl.src.chemecho.get_spectrum import extract_spectrum_data
from ir_visualisation import ir_graph
from ir_spectra_conversion import spectrum_to_dataframe, from_df_to_csv
from musification import spectrum_to_midi

st.title("ChemEcho")
st.write("Convert an IR spectrum into music. Enter a CAS number below.")

cas = st.text_input("CAS number", value="67-64-1", placeholder="e.g. 67-64-1 for acetone")

if st.button("Generate"):
    with st.spinner("Fetching spectrum from NIST..."):
        try:
            compound = nist.get_compound(cas)
            st.subheader(compound.name)

            data = extract_spectrum_data(cas)

            # IR plot
            fig, ax = ir_graph(data, cas)
            st.pyplot(fig)

            # CSV download
            df = spectrum_to_dataframe(data)
            csv_str = from_df_to_csv(df)
            st.download_button(
                label="Download spectrum as CSV",
                data=csv_str,
                file_name=f"{cas}_spectrum.csv",
                mime="text/csv"
            )

            # MIDI generation
            with tempfile.NamedTemporaryFile(suffix=".mid", delete=False) as tmp:
                midi_path = tmp.name

            spectrum_to_midi(cas, midi_path)

            with open(midi_path, "rb") as f:
                st.download_button(
                    label="Download MIDI",
                    data=f,
                    file_name=f"{cas}.mid",
                    mime="audio/midi"
                )
            os.remove(midi_path)

        except Exception as e:
            st.error(f"Something went wrong: {e}")
