from setuptools import find_packages, setup


setup(
    name="ml-stability",
    version="0.1.0",
    author="Janosh Riebesell",
    author_email="janosh.riebesell@gmail.com",
    url="https://github.com/janosh/ml-stability",
    description="Data-driven search for novel stable materials",
    long_description=open("readme.md").read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    keywords=[
        "data-driven materials discovery",
        "finding new stable crystal structures",
        "machine learning",
        "high-throughput",
    ],
    install_requires=[
        "matplotlib",
        "pymatviz",
        "plotly",
        "tqdm",
        "scikit_learn",
    ],
    extras_require={
        "wren": ["aviary"],  # not on PyPI, install manually
        "single_use_deps": [
            "seaborn",
        ],
    },
)
