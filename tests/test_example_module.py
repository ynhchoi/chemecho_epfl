import sys
from pathlib import Path

# Ajoute le chemin du dossier "chemecho_epfl" à sys.path
# Path(__file__).parent = dossier_contenant_test
# Path(__file__).parent.parent = chemecho_epfl
chemin_chemecho_epfl = str(Path(__file__).parent.parent.resolve())
sys.path.append(chemin_chemecho_epfl)

# Maintenant tu peux importer depuis chemecho_epfl_main_code
from chemecho_epfl_main_code.ir_spectra_conversion import from_jdx_to_df
from chemecho_epfl_main_code.ir_visualisation import ir_graph


# Test the function
ir_graph(from_jdx_to_df("7732-18-5-IR.jdx"), "Water")[0].savefig("test.png", bbox_inches="tight", dpi=200)
