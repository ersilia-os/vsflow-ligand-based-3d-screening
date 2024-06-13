from setuptools import setup, find_packages

with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

setup(
    name="vsflow3d",
    version="0.1.0",
    author="Miquel Duran-Frigola, Ersilia Open Source Initiative",
    author_email="miquel@ersilia.io",
    url="https://github.com/ersilia-os/vsflow-ligand-based-3d-screening",
    description="VSFlow based module to obtain 3D shape similarity measures",
    license="GPLv3",
    python_requires=">=3.9",
    install_requires=install_requires,
    extras_require={},
    packages=find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ),
    keywords="3dshape",
    project_urls={"Source Code": "https://github.com/ersilia-os/vsflow-ligand-based-3d-screening/"},
)