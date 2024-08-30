from setuptools import setup, find_packages

setup(
    name='drMD',
    version='1.0.0',
    description='drMD: Molecular Dynamics for Protein Scientists',
    author="Dr Eugene Shrimpton-Phoenix",
    author_email="eshrimpt@ed.ac.uk",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ESPhoenix/drMD",
    packages=find_packages(),
    python_requires='>=3.9',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Linux',
    ],
    install_requires=[
        'pyyaml',
        'argparse',
        'pandas',
        'mdtraj',
        'numpy',
        'matplotlib',
        'weasyprint',
        'scipy',
        'tqdm',
        'pdbutils',
        'scikit-learn',
        'propka',
        'mdanalysis'
    ],
)