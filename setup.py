from setuptools import find_packages, setup

setup(
    name='alpaca',
    packages=find_packages(include=['alpaca_proteomics']),
    version='0.9.98',
    description='Absolute Protein Quantification python library',
    author='Borja Ferrero Bordera',
    author_email='ferrerobob@uni-greifswald.de',
    project_urls={'Author':'https://www.linkedin.com/in/borjaferrero/'},
    license='MIT',
    install_requires= [
            'matplotlib',
            'numpy',
            'pandas',
            'scikit-learn',
            'scipy',
            'seaborn',
            'XlsxWriter']
)
