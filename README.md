![Project Logo](assets/banner.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
Chem Echo
</h1>

<br>


None

## 🔥 Discover Chem Echo, the app that makes molecules sing! 🎶

Chem Echo is a tool to translate IR spectrum into music, making accessible the data for blind or visually impaired person. This package is also intended for curious people who wants to discover a new fun and cool approach to compare molecules. 

To maximize the information contained in the audio file, ChemEcho does not only provide a literal translation but also:
- choice of instrument based on the logP (hence polarity) of the molecule of interest
- tempo based on the molecular weight, a reminder that the lighter the molecule, the faster the vibrations. 

## 👩‍💻 Installation

Create a new environment, you may also give the environment a different name. You can activate it

```
conda create -n chemecho python=3.10 
conda activate chemecho
```
Then, install the package using one of the following way. Whether using pip:
```
(chemecho) $ pip install chemecho
```
Or by cloning the GitHub repository:
```
pip install git+github.com:yhchoi/chemecho.git 
```

If somehow you miss one or several packages needed for the app after the installation (they should be automatically installed), here is the list of the requirements and how to install them (if needed). 
``` 
pip install matplotlib streamlit jcamp musicpy nistchempy io os pandas requests tempfile base6 scipy rdkit 
```

An additional module (daw) needs to be installed for the musicpy package. To do so, visit the link: https://musicpy.readthedocs.io/en/latest/Musicpy%20daw%20module/. Read and follow the instruction in the "Preparation before importing" paragraph. 

ChemEcho possesses its own streamlit app. However, if you want to test some functions independantly from the app, you may need a jupyter notebook. If you need jupyter lab, install it 

```
(chemecho) $ pip install jupyterlab
```

## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:yhchoi/chemecho`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:yhchoi/chemecho.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(chemecho) $ pip install -e ".[test,doc]"
```

As for the classic installation, you will need to install the daw module from musicpy package. 

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```


```python
from mypackage import main_func

# One line to rule them all
result = main_func(data)
```

This usage example shows how to quickly leverage the package's main functionality with just one line of code (or a few lines of code). 
After importing the `main_func` (to be renamed by you), you simply pass in your `data` and get the `result` (this is just an example, your package might have other inputs and outputs). 
Short and sweet, but the real power lies in the detailed documentation.



