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
    view.setBackgroundColor("white")
    view.zoomTo()
    return view._make_html()


# ── Page config (must be first Streamlit call) ────────────────────────────────

st.set_page_config(
    page_title="ChemEcho",
    page_icon="🎶",
    layout="wide",
    menu_items={
        "About": "ChemEcho — sonification of IR spectra for accessible chemistry.",
    },
)


# ── Sidebar ───────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("# About")
    st.markdown(
        "ChemEcho turns an IR spectrum into music so the structure and "
        "functional groups of a molecule can be perceived through sound."
    )
    st.markdown("---")
    st.markdown("**Sound layers**")
    st.markdown(
        "- 🥁 Carbon prelude — one drum hit per carbon atom\n"
        "- 🎻 Cello melody — traces the IR transmittance contour\n"
        "- 🎺 Functional-group accents — distinctive instrument per FG\n"
        "- 🎼 Tempo — scales with molecular weight (lighter = faster)"
    )


# ── UI ────────────────────────────────────────────────────────────────────────

st.title("Chem Echo 🎶")
st.write(
    "See a molecule's structure, IR spectrum, and hear it as music. "
    "Designed to make IR spectroscopy accessible through sound."
)
st.markdown(
    "**Enter *any one* of the following — just one is enough:**\n"
    "- 🏷️ **Name** — e.g. `caffeine`, `acetone`, `aspirin`\n"
    "- 🔢 **CAS number** — e.g. `67-64-1` (acetone)\n"
    "- 🧪 **Molecular formula** — e.g. `C6H6` (benzene)\n"
    "- 🧬 **SMILES** — e.g. `CC(=O)C` (acetone)"
)

with st.expander("ℹ️ How to listen — what each sound means"):
    st.markdown(
        """
The music is built from your molecule's infrared (IR) spectrum, in four layers
that play together:

**1. Carbon prelude** (the very beginning)
A short drum sequence — one hit per carbon atom — counts the molecule's
carbon skeleton before the spectrum starts.

**2. Cello main melody**
A cello traces the IR spectrum shape. Deeper absorption peaks (lower
transmittance) produce higher pitches. The melody plays from high wavenumber
to low — the conventional left-to-right direction of an IR plot.

**3. Wavenumber ruler (background drums)**
A snare drum marks every 100 cm⁻¹, a kick drum every 500 cm⁻¹. These act
as your "axis ticks" so you can orient where you are in the spectrum.

**4. Functional-group accents**
When a known IR absorption band is both predicted from the structure
(SMARTS pattern) **and** confirmed by a real peak in the spectrum, a
distinctive instrument plays a signature pitch at that moment:

| Functional group | Instrument |
|---|---|
| O–H (alcohol) | piccolo |
| O–H (carboxylic acid) | clarinet |
| N–H | oboe |
| C–H (saturated) | nylon guitar |
| C–H (aromatic) | steel guitar |
| C≡N (nitrile) | music box |
| C≡C (alkyne) | vibraphone |
| C=O (carbonyl) | orchestra hit |
| C=C (alkene) | violin |
| aromatic ring | flute |
| C–O | French horn |
| N=O (nitro) | trumpet |

Louder accent = deeper (more prominent) IR peak.

**Tempo encodes molecular weight.** Lighter molecules play faster, heavier
slower — derived from kinetic theory (Maxwell-Boltzmann mean speed),
bounded to 60–120 BPM so the music always stays intelligible.
        """
    )

query = st.text_input(
    "Molecule",
    placeholder="e.g. caffeine  |  67-64-1  |  C6H6  |  CC(=O)C",
    help="Accepts a compound name, CAS number, molecular formula, or SMILES string.",
)

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
                            st.image(
                                img,
                                caption=f"2D skeletal structure of {compound.name}",
                            )
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
                wn_min = int(min(extracted_data[0]))
                wn_max = int(max(extracted_data[0]))
                st.caption(
                    f"IR transmittance spectrum of {compound.name}, "
                    f"covering {wn_min} to {wn_max} cm⁻¹. "
                    f"{len(extracted_data[0])} data points."
                )
                buffer = io.StringIO()
                fig.savefig(buffer, format="svg")
                svg = buffer.getvalue()

                st.download_button(
                    label="Download spectrum as SVG",
                    data=svg,
                    file_name=f"{compound.name}_spectrum.svg",
                    mime="image/svg+xml")

            if not smiles:
                raise ValueError(
                    "SMILES not available — functional-group detection requires "
                    "a valid SMILES from PubChem."
                )

            # MIDI generation
            midi_path = None
            legend = None
            try:
                midi_path, legend = molecular_music_fg(extracted_data, compound, smiles)
                with open(midi_path, "rb") as f:
                    midi_bytes = f.read()
            finally:
                if midi_path and os.path.exists(midi_path):
                    os.remove(midi_path)

            st.markdown("**Listen**")

            # Context-aware "what you're about to hear" — molecule-specific.
            # Helps blind / VI users and first-time listeners orient before playback.
            preview_parts = []
            if legend['carbon_count'] > 0:
                preview_parts.append(
                    f"{legend['carbon_count']} drum hit"
                    f"{'s' if legend['carbon_count'] != 1 else ''} for carbon count"
                )
            preview_parts.append(f"cello melody at {legend['bpm']} BPM")
            n_fgs = len(legend['fgs'])
            if n_fgs > 0:
                preview_parts.append(
                    f"{n_fgs} functional-group accent"
                    f"{'s' if n_fgs != 1 else ''}"
                )
            st.info("🎧 You'll hear: " + " → ".join(preview_parts) + ".")

            # Web-friendly GM SoundFont (Magenta format). html-midi-player only
            # supports Magenta-format banks, so this is currently the best
            # publicly-hosted general-MIDI option. The same MIDI file will
            # sound different in desktop DAWs (GarageBand, Logic, etc.) because
            # those use their own, much larger instrument libraries.
            _SOUNDFONT_URL = (
                "https://storage.googleapis.com/magentadata/js/soundfonts/sgm_plus"
            )
            midi_b64 = base64.b64encode(midi_bytes).decode()
            components.html(f"""
                <script src="https://cdn.jsdelivr.net/npm/tone@14/build/Tone.js"></script>
                <script src="https://cdn.jsdelivr.net/npm/@magenta/music@1.23.1/es6/core.js"></script>
                <script src="https://cdn.jsdelivr.net/npm/html-midi-player@1.5.0"></script>
                <midi-player
                    src="data:audio/midi;base64,{midi_b64}"
                    sound-font="{_SOUNDFONT_URL}"
                    style="width:100%;margin-top:4px"
                    aria-label="MIDI player for {compound.name} spectrum sonification">
                </midi-player>
            """, height=150)

            st.download_button(
                label="Download MIDI",
                data=midi_bytes,
                file_name=f"{compound.name}.mid",
                mime="audio/midi"
            )

            # FG legend (v2 only) — presented BOTH as a textual list (for screen
            # readers) AND as a table (for sighted users).
            if legend:
                st.markdown("**Detected functional groups**")
                st.markdown(
                    f"Molecular weight: **{compound.mol_weight:.2f} g/mol**  ·  "
                    f"Carbon count: **{legend['carbon_count']}**  ·  "
                    f"Tempo: **{legend['bpm']} BPM**"
                )

                if legend['fgs']:
                    st.markdown(
                        f"Found **{len(legend['fgs'])}** functional group"
                        f"{'s' if len(legend['fgs']) != 1 else ''}:"
                    )
                    for fg, info in legend['fgs'].items():
                        lo, hi = info['region']
                        n = info['n_peaks']
                        st.markdown(
                            f"- **{fg}** — {n} peak{'s' if n != 1 else ''} "
                            f"in {lo}–{hi} cm⁻¹"
                        )

                    rows = [
                        {
                            "Functional group": fg,
                            "Region (cm⁻¹)": f"{info['region'][0]}–{info['region'][1]}",
                            "Peaks detected": info['n_peaks'],
                        }
                        for fg, info in legend['fgs'].items()
                    ]
                    st.table(pd.DataFrame(rows))
                else:
                    st.markdown(
                        "No functional groups confirmed in this spectrum "
                        "(no IR peaks fell in any catalog region)."
                    )

        except Exception as e:
            st.error(f"Something went wrong: {e}")
