
from setuptools import setup, find_packages

setup(
    name="drakkar",
    version="1.0.0",
    author="Antton Alberdi",
    author_email="antton.alberdi@sund.ku.dk",
    description="Metagenomics pipeline optimised for Mjolnir",
    packages=find_packages(),
    py_modules=["drakkar"],
    install_requires=[
        "numpy",
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "drakkar=drakkar.cli:main",
        ],
    },
    python_requires=">=3.6",
)
