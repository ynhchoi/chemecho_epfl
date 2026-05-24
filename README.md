![Project Logo](assets/banner.png)

<h1 align="center">ChemEcho 🎶</h1>



## 🔥 Discover ChemEcho, the app that makes molecules sing!

ChemEcho is a tool that translates a molecule's infrared (IR) spectrum into
music, making spectroscopy data accessible to blind and visually impaired
users. It is also intended for curious people who want a new, fun way to
explore and compare molecules by ear.

> ⚠️ **Internet access is required.** ChemEcho looks up molecule
> structures on **PubChem** and fetches IR spectra from the **NIST WebBook**
> at run time. The package cannot be used offline.

To maximize the information contained in the audio file, ChemEcho does not
only provide a literal translation of the spectrum, but layers four cues:

| Layer | Encodes | How |
|---|---|---|
| Main melody | IR transmittance curve | Pitch ∝ (100 − transmittance); rendered on a Fantasia synth pad |
| Functional-group notes | Detected groups (C=O, O–H, ...) | One instrument per group, fired when an absorption peak falls in its region |
| Prelude drum hits | Carbon count | One Synth-Drum hit per carbon, played before the melody |
| Tempo (BPM) | Molecular weight | Maxwell–Boltzmann mapping (v ∝ 1/√M): lighter ⇒ faster |

Molecules can be entered in any of four formats — **name, CAS number,
formula, or SMILES** — and ChemEcho auto-detects which.

## 👩‍💻 Installation

ChemEcho supports **Python 3.10** on macOS, Linux, and
Windows

Create a new environment (you may give it any name you like) and activate
it. The same `conda` commands work on macOS, Linux, and Windows — use a
regular terminal on macOS/Linux, and **Anaconda Prompt** on Windows:

```
conda create -n chemecho python=3.10
conda activate chemecho
```

Then install the package from the GitHub repository:

```
(chemecho) $ pip install git+https://github.com/ynhchoi/ChemEcho.git
```

**Verify that the installation worked:**

```
(chemecho) $ python -c "import chemecho; print('ChemEcho', chemecho.__version__, 'installed OK')"
```

If somehow one or several packages required by the app are missing after
installation (they should be installed automatically), here is the manual
fallback:

```
pip install rdkit matplotlib "nistchempy>=2.0.0" pandas streamlit musicpy jcamp requests "numpy>=1.25,<2.0" scipy py3Dmol Pillow
```

### Audio rendering — `musicpy.daw`

An additional module (`daw`) is needed by the `musicpy` package for audio
rendering. Follow the instructions in the "Preparation before importing"
paragraph here:
<https://musicpy.readthedocs.io/en/latest/Musicpy%20daw%20module/>


### Launch Streamlit app

ChemEcho ships its own Streamlit app. 

Run from the project root:

```

(chemecho) $ streamlit run src/chemecho/streamlit_app.py
```


If you want to experiment with the functions outside the app, install
JupyterLab and open the example notebook in `notebooks/report.ipynb`:

```
(chemecho) $ pip install jupyterlab
(chemecho) $ jupyter lab
```

## 🛠️ Development installation

For working on the ChemEcho codebase itself.

### 1. Get the source

If you want to contribute back to this repository:

```
git clone https://github.com/ynhchoi/ChemEcho.git
cd ChemEcho
```

If you want to fork the project under your own GitHub account, create an
empty repository on `https://github.com/<your-username>/ChemEcho`
first, then initialize and push (works on macOS, Linux, and Windows; use
**Anaconda Prompt**, **Git Bash**, or any modern terminal on Windows):

```
git init
git add .
git commit -m "Initial commit"
git branch -M main
git remote add origin git@github.com:<your-username>/ChemEcho.git
git push -u origin main
```

### 2. Install in editable mode

Install the package together with the test and documentation extras:

```
(chemecho) $ pip install -e ".[test,doc]"
```

As with the regular installation, you still need to install the `daw`
module from the `musicpy` package (see § *Audio rendering — `musicpy.daw`*
above).

### 3. Run tests and coverage

```
(chemecho) $ pip install tox
(chemecho) $ tox
```

## 🧯 Troubleshooting

### `pip install chemecho` fails with "No matching distribution found"
ChemEcho is not published on PyPI, so plain `pip install chemecho` will not work.
Install directly from GitHub using the `pip install git+...` command
in the Installation section above.

### `ValueError: PubChem could not find '<your input>'`
PubChem did not recognize your input. Double-check the spelling of the
compound name, the CAS number format (`12345-67-8`), or try a different
identifier (e.g. switch from name to SMILES).

### `ValueError: No IR data found on NIST for '<molecule>'`
NIST WebBook does not have an IR spectrum for this molecule. Try a closely
related compound — not every molecule in the database has IR data.

### `requests.exceptions.ConnectionError` or `Timeout`
Network problem reaching PubChem or NIST. Check your internet connection
(remember: ChemEcho cannot work offline) and/or retry after a few seconds.

### The Streamlit app shows **"No peaks detected"**
The spectrum may be very flat or noisy. ChemEcho filters out peaks below
5 % of the transmittance range by design.

## 👥 Authors

- Eulalie Schwendenmann · `eulalie.schwendenmann@epfl.ch`
- Zoé Chartier · `zoe.chartier@epfl.ch`
- Yunhee Choi · `yunhee.choi@epfl.ch`

EPFL · CH-200 *Practical programming in chemistry*

## 📄 License

Released under the MIT License — see [`LICENSE`](LICENSE).
