from setuptools import setup

from matbench_discovery import ROOT, URLs

setup(
    name="matbench-discovery",
    version="0.1.0",
    author="Janosh Riebesell, Rhys Goodall",
    author_email="janosh@lbl.gov",
    url="https://github.com/janosh/matbench-discovery",
    description="A machine learning benchmark that simulates high-throughput screening "
    "for new materials and ranks energy models by their ability to increase the hit "
    "rate of stable crystals",
    long_description=open("readme.md").read(),
    long_description_content_type="text/markdown",
    packages=["matbench_discovery"],
    python_requires=">=3.9",
    package_data={
        "matbench_discovery": [f"{ROOT}/data/mp/*.json"],
    },
    keywords=[
        "data-driven materials discovery",
        "crystal stability",
        "machine learning",
        "materials space",
        "high-throughput search",
        "energy above convex hull",
        "formation energy",
    ],
    install_requires=[
        "matplotlib",
        "pymatgen",
        "numpy",
        "pandas",
        "scikit-learn",
        "scipy",
        "plotly",
        "tqdm",
        "wandb",
    ],
    extras_require={
        "test": ["pytest", "pytest-cov", "pytest-markdown-docs"],
        "running-models": ["aviary", "m3gnet", "maml", "megnet"],
    },
    project_urls=URLs,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
