![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)
![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)

### Pairs-based analysis toolkit

### Prepare environment
1. make sure that You have python >= 3.10 installed as well as R >= 4.0.0
2. install [poetry](https://python-poetry.org/) dependency manager

       pip install poetry

3. clone repository

       git clone https://github.com/ClinicalEpigeneticsLaboratory/PairsBasedAnalysis.git

4. open project directory and install required dependencies

       poetry install 
       Rscript r_install.R

### External tools necessary to run scripts:

1. HOMER -> [link to manual](http://homer.ucsd.edu/homer/)
2. LOLA DATABASE -> [instruction to build database](https://databio.org/regiondb)

### Examples of usage are placed in:
1. Step 1: Identification of CpG pairs - example_1.ipynb
2. Step 2: Analysis of identified CpG pairs - example_2.ipynb


To run these `notebooks` start jupyterlab session:
    
    poetry run python -m jupyterlab