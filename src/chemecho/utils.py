import re
import time
import nistchempy as nist
import requests
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol
from IPython.display import HTML

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
    return HTML(view._make_html())