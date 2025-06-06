from setuptools import setup, find_packages

setup(
    name="drakkar",
    version="1.0.0",
    author="Antton Alberdi",
    author_email="antton.alberdi@sund.ku.dk",
    description="Metagenomics pipeline optimized for Mjolnir",
    packages=find_packages(),  # Automatically discovers drakkar/
    include_package_data=True,  # Ensures non-Python files are included
    install_requires=[
        "numpy",
        "pandas",
        "argparse",
        "PyYAML",
        "plotly",
        "biopython"
    ],
    entry_points={
        "console_scripts": [
            "drakkar=drakkar.cli:main",
        ],
    },
    python_requires=">=3.6",
)
