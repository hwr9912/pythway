from setuptools import setup, find_packages

setup(
    name="pythway",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "requests",
        "lxml",
    ],
    entry_points={
        "console_scripts": [
            "kegg_gmt=my_kegg_gmt.kegg_gmt:main",
        ],
    },
)
