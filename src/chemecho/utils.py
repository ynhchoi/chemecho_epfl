import re
import time
import requests
import py3Dmol
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import nistchempy as nist


_CAS_RE = re.compile(r'^\d{2,7}-\d{2}-\d$')
_FORMULA_RE = re.compile(r'^([A-Z][a-z]?\d*)+$')
_SMILES_CHARS = set('=#()[]@/\\+%')
_PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"


def _detect_input_type(query: str) -> str:
    """
    Detects the type of chemical identifier in input

    Args:
        query (str): user input string (CAS number, SMILES, molecular formula, or name)

    Return:
        (str): one of 'cas', 'smiles', 'formula', or 'name'
    """
    if _CAS_RE.match(query):
        return 'cas'
    if any(c in query for c in _SMILES_CHARS):
        return 'smiles'
    if _FORMULA_RE.match(query):
        return 'formula'
    return 'name'


def _pubchem_props(url: str) -> dict:
    """
    Fetches compound properties from a PubChem REST API URL.

    Args:
        url (str): full PubChem property URL to query

    Return:
        (dict): first Properties entry from the PubChem response, or empty dict if
                the request fails or returns a non-200 status code
    """
    resp = requests.get(url, timeout=10)
    if resp.status_code != 200:
        return {}
    return resp.json().get('PropertyTable', {}).get('Properties', [{}])[0]


def _formula_props(formula: str) -> dict:
    """
    Fetches compound properties from PubChem using molecular formula
    Uses asynchronous ListKey polling bc formula searches are not immediate.

    Args:
        formula (str): molecular formula of the compound (e.g. 'C6H12O6')

    Return:
        (dict): first Properties entry from the PubChem response, or empty dict if
                the request fails, returns a non-202 status, or times out after ~15 seconds
    """
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
    """
    Extracts a SMILES string from a PubChem properties dictionary.
    Tries multiple key name variants in order of preference.

    Args:
        props (dict): PubChem properties dictionary

    Return:
        (str | None): SMILES string if found, None if no SMILES key is present
    """
    return (props.get('IsomericSMILES') or props.get('CanonicalSMILES')
            or props.get('SMILES') or props.get('ConnectivitySMILES'))


def resolve_molecule(query: str) -> dict:
    """
    Resolves a chemical identifier to its SMILES, name, and CAS number via PubChem.
    Accepts CAS number, SMILES, molecular formula, and name

    Args:
        query (str): chemical identifier (CAS number, SMILES, molecular formula, or name)

    Return:
        (dict): dictionary with keys:
                - 'smiles' (str): SMILES string of the compound
                - 'name' (str): IUPAC name of the compound
                - 'cas' (str | None): CAS number if found in PubChem synonyms, else None

    Raises:
        ValueError: if PubChem cannot find the compound
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
    """
    Fetches a NIST compound object with IR data, trying CAS number first
    then falling back to a name search if the CAS fails or is not provided.

    Args:
        cas (str | None): CAS number of the compound, or None to skip directly to name search
        name (str): name of the compound, used as fallback if CAS lookup fails

    Return:
        (tuple): (cas (str), compound (NistCompound)) where cas is the CAS number
                 used to retrieve the compound

    Raises:
        ValueError: if no IR data is found on NIST for the given name
    """
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
    """
    Generates a 2D structural image of a molecule from its SMILES

    Args:
        smiles (str): SMILES of compound

    Return:
        (PIL.Image | None): 400x300 PIL image of the molecule, or None if the
                            SMILES string is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Draw.MolToImage(mol, size=(400, 300))


def molecule_3d_html(smiles: str) -> str | None:
    """
    Generates an interactive 3D molecular viewer as an HTML string using py3Dmol.
    3D coordinates are computed with RDKit using MMFF force field optimization

    Args:
        smiles (str): SMILES string of the compound

    Return:
        (str | None): HTML string containing the interactive 3D viewer, or None if
                      the SMILES is invalid or 3D coordinate generation fails
    """
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