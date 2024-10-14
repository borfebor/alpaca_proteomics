from setuptools import find_packages, setup

setup(
    name='alpaca proteomics',
    packages=find_packages(include=['alpaca_proteomics']),
    version='0.9.98.240922',
    description='Absolute Protein Quantification python library',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Borja Ferrero Bordera',
    author_email='b.ferrero@med.uni-muenchen.de',
    project_urls={'Author':'https://www.linkedin.com/in/borjaferrero/'},
    license='MIT',
    install_requires= [
            'matplotlib',
            'numpy',
            'pandas',
            'scikit-learn',
            'scipy',
            'seaborn',
            'thefuzz',
            'XlsxWriter']
)
