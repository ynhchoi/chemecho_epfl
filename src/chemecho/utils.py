import nistchempy as nist
import requests
from rdkit import Chem
from rdkit.Chem import Draw

def name_to_cas(name: str) -> str:
    """Search NIST by compound name and return CAS number.

    Args:
        name (str): name of the compound of interest

    Returns:
        (str): corresponding CAS number    
    """

    results = nist.run_search(name, 'name', cIR=True)
    if not results.compounds:
        raise ValueError(f"No compound with IR data found for '{name}' on NIST.")
    
    compound = results.compounds[0]
    cas = compound.cas_rn
    return cas


def get_smiles(cas: str) -> str | None:
    """Fetch canonical SMILES from PubChem by compound name.
    
    Args:
        cas (str): CAS number of molecule of interest

    Returns:
        (str): corresponding SMILES
    """
    
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{cas}/property/IsomericSMILES/JSON")
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            props = resp.json()['PropertyTable']['Properties'][0]
            return props.get('IsomericSMILES') or props.get('SMILES')
    except Exception:
        pass
    return None


def draw_molecule(smiles: str):
    """Return a PIL image of the molecule from its SMILES string.
    
    Args:
        smiles (str): SMILES of the molecule of interest
        
    Returns:
        (PIL.Image.Image): corresponding PIL image 
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Draw.MolToImage(mol, size=(400, 300))
