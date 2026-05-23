import re
import time
import streamlit as st
import streamlit.components.v1 as components
import nistchempy as nist
import pandas as pd
import requests
import os
import io
import base64
import py3Dmol
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

from get_spectrum import extract_spectrum_data, ir_graph, from_df_to_csv
from musification import molecular_music
from musification_fg import molecular_music_fg


_CAS_RE = re.compile(r'^\d{2,7}-\d{2}-\d$')
_FORMULA_RE = re.compile(r'^([A-Z][a-z]?\d*)+$')
_SMILES_CHARS = set('=#()[]@/\\+%')
_PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"


def _detect_input_type(query: str) -> str:
    if _CAS_RE.match(query):
        return 'cas'
    if any(c in query for c in _SMILES_CHARS):
        return 'smiles'
    if _FORMULA_RE.match(query):
        return 'formula'
    return 'name'


def _pubchem_props(url: str) -> dict:
    """GET a PubChem property URL and return the first Properties entry."""
    resp = requests.get(url, timeout=10)
    if resp.status_code != 200:
        return {}
    return resp.json().get('PropertyTable', {}).get('Properties', [{}])[0]


def _formula_props(formula: str) -> dict:
    """Async formula search via PubChem ListKey polling (max ~15 s)."""
    r = requests.get(
        f"{_PUBCHEM}/formula/{formula}/JSON?MaxRecords=1", timeout=10
    )
    if r.status_code != 202:
        return {}
    listkey = r.json()['Waiting']['ListKey']
    poll = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{listkey}/property/CanonicalSMILES,IUPACName/JSON"
    for _ in range(8):
        time.sleep(2)
        r = requests.get(poll, timeout=10)
        if r.status_code == 200:
            return r.json().get('PropertyTable', {}).get('Properties', [{}])[0]
    return {}


def _smiles_from(props: dict) -> str | None:
    """Extract SMILES from PubChem props regardless of key name variant."""
    return (props.get('IsomericSMILES') or props.get('CanonicalSMILES')
            or props.get('SMILES') or props.get('ConnectivitySMILES'))


def resolve_molecule(query: str) -> dict:
    """Resolve name / CAS / formula / SMILES via PubChem.

    Returns dict with keys: smiles, name, cas (cas may be None).
    Raises ValueError if PubChem cannot find the compound.
    """
    q = query.strip()
    kind = _detect_input_type(q)
    _props = "CanonicalSMILES,IsomericSMILES,IUPACName"

    if kind == 'smiles':
        props = _pubchem_props(f"{_PUBCHEM}/smiles/{requests.utils.quote(q)}/property/{_props}/JSON")
    elif kind == 'formula':
        props = _pubchem_props(f"{_PUBCHEM}/name/{q}/property/{_props}/JSON")
        if not props:
            props = _formula_props(q)
    else:
        props = _pubchem_props(f"{_PUBCHEM}/name/{requests.utils.quote(q)}/property/{_props}/JSON")

    if not props:
        raise ValueError(f"PubChem could not find '{q}' (detected as: {kind}).")

    cid = props['CID']
    smiles = _smiles_from(props)
    name = props.get('IUPACName') or q

    cas = None
    syns_resp = requests.get(f"{_PUBCHEM}/cid/{cid}/synonyms/JSON", timeout=10)
    if syns_resp.status_code == 200:
        synonyms = syns_resp.json()['InformationList']['Information'][0].get('Synonym', [])
        for syn in synonyms:
            if _CAS_RE.match(syn):
                cas = syn
                break

    return {'smiles': smiles, 'name': name, 'cas': cas}


def nist_compound_from(cas: str | None, name: str):
    """Fetch NIST compound with IR data, trying CAS first then name search."""
    if cas:
        try:
            return cas, nist.get_compound(cas)
        except Exception:
            pass

    results = nist.run_search(name, 'name', cIR=True)
    if not results.compounds:
        raise ValueError(f"No IR data found on NIST for '{name}'.")
    nist_id = results.compounds[0].ID
    numeric = nist_id.lstrip('C').split('&')[0]
    cas = f"{numeric[:-3]}-{numeric[-3:-1]}-{numeric[-1]}"
    return cas, nist.get_compound(cas)


def draw_molecule(smiles: str):
    """Return a PIL image of the molecule from its SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Draw.MolToImage(mol, size=(400, 300))


def molecule_3d_html(smiles: str) -> str | None:
    """Generate 3D coordinates and return an HTML string for py3Dmol viewer."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, randomSeed=42)
    if result != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol)
    molblock = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width="100%", height=380)
    view.addModel(molblock, "mol")
    view.setStyle({"stick": {"colorscheme": "Jmol"}, "sphere": {"scale": 0.25, "colorscheme": "Jmol"}})
    view.setBackgroundColor("#0e1117")
    view.zoomTo()
    view.spin(True)
    return view._make_html()


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
                st.markdown("**Detected functional groups**")
                st.caption(f"BPM: {legend['bpm']}  ·  Carbon count: {legend['carbon_count']}")
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
