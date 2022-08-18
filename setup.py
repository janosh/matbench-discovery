from setuptools import find_packages, setup


setup(
    name="mb-discovery",
    version="0.1.0",
    author="Janosh Riebesell, Rhys Goodall",
    author_email="janosh.riebesell@gmail.com",
    url="https://github.com/janosh/matbench-discovery",
    description="Data-driven search for novel stable materials",
    long_description=open("readme.md").read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    keywords=[
        "data-driven materials discovery",
        "finding new stable crystal structures",
        "machine learning",
        "high-throughput",
        "energy above convex hull",
    ],
    install_requires=[
        "matplotlib",
        "pymatgen",
        "numpy",
        "pandas",
        "scipy",
        "plotly",
        "tqdm",
    ],
    extras_require={
        "running-models": ["wandb", "m3gnet"],
    },
)
