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

ChemEcho supports **Python 3.10, 3.11, and 3.12** on macOS, Linux, and
Windows. (Python 3.10 recommended)

Create a new environment (you may give it any name you like) and activate
it. The same `conda` commands work on macOS, Linux, and Windows — use a
regular terminal on macOS/Linux, and **Anaconda Prompt** on Windows:

```
conda create -n chemecho python=3.10
conda activate chemecho
```

Then install the package. Either from PyPI:

```
(chemecho) $ pip install chemecho
```

Or directly from the GitHub repository:

```
(chemecho) $ pip install git+https://github.com/ynhchoi/chemecho_epfl.git
```

**Verify that the installation worked:**

```
(chemecho) $ python -c "import chemecho; print('ChemEcho', chemecho.__version__, 'installed OK')"
```

If somehow one or several packages required by the app are missing after
installation (they should be installed automatically), here is the manual
fallback:

```
pip install rdkit matplotlib nistchempy==1.0.5 pandas streamlit musicpy jcamp requests "numpy>=1.25,<2.0" scipy py3Dmol
```

### Audio rendering — `musicpy.daw`

An additional module (`daw`) is needed by the `musicpy` package for audio
rendering. Follow the instructions in the "Preparation before importing"
paragraph here:
<https://musicpy.readthedocs.io/en/latest/Musicpy%20daw%20module/>

`musicpy.daw` relies on **FluidSynth**. Install it for your OS, then make
sure Python can find the shared library:

**macOS (Homebrew):**

```
brew install fluid-synth
export DYLD_LIBRARY_PATH=/opt/homebrew/lib
```

(Add the `export` line to `~/.zshrc` or `~/.bash_profile` to make it permanent.)

**Linux (Ubuntu/Debian):**

```
sudo apt-get update
sudo apt-get install -y fluidsynth libfluidsynth-dev
```

If `musicpy.daw` still cannot find the library, point the loader at it:

```
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
```

(Add it to `~/.bashrc` to make it permanent.)

**Windows:**

Download a FluidSynth release from
<https://github.com/FluidSynth/fluidsynth/releases>, unzip it, and add the
folder containing `libfluidsynth-*.dll` to your `PATH`:

- *PowerShell (current session):* `$env:PATH = "C:\path\to\fluidsynth\bin;" + $env:PATH`
- *Permanent:* System Properties → Environment Variables → edit `Path`.

### Streamlit app

ChemEcho ships its own Streamlit app. How to launch it depends on how you
installed:

**If you cloned the repository** (run from the project root):

```
(chemecho) $ streamlit run src/chemecho/streamlit_app.py
```

**If you installed from PyPI** (the app file lives inside `site-packages`,
so let Python locate it for you):

```
(chemecho) $ streamlit run "$(python -c 'import chemecho, os; print(os.path.join(os.path.dirname(chemecho.__file__), "streamlit_app.py"))')"
```

On Windows PowerShell, the equivalent one-liner is:

```
(chemecho) $ streamlit run (python -c "import chemecho, os; print(os.path.join(os.path.dirname(chemecho.__file__), 'streamlit_app.py'))")
```

If you want to experiment with the functions outside the app, install
JupyterLab and open the example notebook in `notebooks/report.ipynb`:

```
(chemecho) $ pip install jupyterlab
```

## 🛠️ Development installation

For working on the ChemEcho codebase itself.

### 1. Get the source

If you want to contribute back to this repository:

```
git clone https://github.com/ynhchoi/chemecho_epfl.git
cd chemecho_epfl
```

If you want to fork the project under your own GitHub account, create an
empty repository on `https://github.com/<your-username>/chemecho_epfl`
first, then initialize and push (works on macOS, Linux, and Windows; use
**Anaconda Prompt**, **Git Bash**, or any modern terminal on Windows):

```
git init
git add .
git commit -m "Initial commit"
git branch -M main
git remote add origin git@github.com:<your-username>/chemecho_epfl.git
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
Update `pip` (`python -m pip install --upgrade pip`) and make sure your
environment uses **Python 3.10, 3.11, or 3.12**. If the problem persists,
install directly from GitHub using the `pip install git+...` command in § 1.

### `ImportError: ... libfluidsynth ...` when running the app
FluidSynth is not on Python's library search path. Re-read the
**Audio rendering — `musicpy.daw`** section above and apply the
platform-specific fix for your OS.

### `ValueError: PubChem could not find '<your input>'`
PubChem did not recognize your input. Double-check the spelling of the
compound name, the CAS number format (`12345-67-8`), or try a different
identifier (e.g. switch from name to SMILES).

### `ValueError: No IR data found on NIST for '<molecule>'`
NIST WebBook does not have an IR spectrum for this molecule. Try a closely
related compound — not every molecule in the database has IR data.

### `requests.exceptions.ConnectionError` or `Timeout`
Network problem reaching PubChem or NIST. Check your internet connection
(remember: ChemEcho cannot work offline) and retry after a few seconds.

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
