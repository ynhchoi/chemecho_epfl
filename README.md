![Project Logo](assets/banner.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
chemecho
</h1>

<br>


None

## 🔥 Discover ChemEcho, the app that makes molecules sing!

ChemEcho is a tool to translate IR spectrum into music, making accessible the data for blind or visually impaired person. This package is also intended for curious people who wants to discover a new fun and cool approach to compare molecules. 

## 👩‍💻 Installation

Create a new environment, you may also give the environment a different name. 

```
conda create -n chemecho python=3.10 
```

```
conda activate chemecho
(conda_env) $ pip install .
```

If you need jupyter lab, install it 

```
(chemecho) $ pip install jupyterlab
```
If somehow you miss one or several packages needed for the app after the installation, here is the list of the requirements: matplotlib, streamlit, jcamp, musicpy, nistchempy

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



