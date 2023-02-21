from setuptools import setup

from matbench_discovery import ROOT, URLs, pkg

setup(
    name="matbench-discovery",
    version="0.1.0",
    author="Janosh Riebesell, Rhys Goodall",
    author_email="janosh@lbl.gov",
    url="https://github.com/janosh/matbench-discovery",
    description="A benchmark for machine learning energy models on inorganic crystal "
    "stability prediction from unrelaxed structures",
    long_description=open("readme.md").read(),  # noqa: SIM115
    long_description_content_type="text/markdown",
    packages=["matbench_discovery"],
    python_requires=">=3.9",
    package_data={
        "matbench_discovery": [f"{ROOT}/data/mp/*.json"],
    },
    keywords=pkg["keywords"],
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
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
