# Getting Started

Welcome to the documentation for alpaca_proteomics 🎉

This guide will help you install the package and set up your environment.

---

## 1. Install the package

```py 
pip install alpaca_proteomics
```

## 2. Import the package


```py
from alpaca_proteomics import alpaca 
```


## 🧩 Requirements

Before installing, make sure you have the following:

- Python 3.8+
- pip
- (Optional) virtualenv

You can check your Python version with:

	bash
	python --version

	git clone https://github.com/borfebor/alpaca_proteomics.git
	cd alpaca_proteomics
	pip install -e .

	python -c "import alpaca_proteomics; print(alpaca_proteomics.__version__)"

The following packages are required for the correct execution of alpaca_proteomics:

- matplotlib
- numpy
- pandas
- scikit-learn
- scipy
- seaborn
- thefuzz
- XlsxWriter
