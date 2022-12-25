from setuptools import find_packages, setup

setup(
    name="matbench-discovery",
    version="0.1.0",
    author="Janosh Riebesell, Rhys Goodall",
    author_email="janosh.riebesell@gmail.com",
    url="https://github.com/janosh/matbench-discovery",
    description="A machine learning benchmark for energy models that emulates high-throughput "
    "materials screening to test which models most accelerate the discovery of new stable "
    "crystals.",
    long_description=open("readme.md").read(),
    long_description_content_type="text/markdown",
    packages=find_packages(exclude=["tests", "tests.*"]),
    python_requires=">=3.9",
    keywords=[
        "data-driven materials discovery",
        "crystal stability",
        "machine learning",
        "high-throughput search",
        "energy above convex hull",
        "formation energy",
    ],
    install_requires=[
        "matplotlib",
        "pymatgen",
        "numpy",
        "pandas",
        "scipy",
        "plotly",
        "tqdm",
        "wandb",
    ],
    extras_require={
        "test": ["pytest", "pytest-cov", "pytest-markdown-docs"],
        "running-models": ["m3gnet", "aviary", "maml", "megnet", "m3gnet-dgl"],
    },
    project_urls={
        "Docs": "https://matbench-discovery.janosh.dev",
        "Package": "https://pypi.org/project/matbench-discovery",
    },
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
