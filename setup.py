import setuptools

with open("readme.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ysospy-noahbraunf",
    version="0.0.2",
    author="Noah Braunfeld",
    author_email="noah@braunfeld.dev",
    description="Useful helper functions for HRL",
    long_description=long_description,
    url="https://github.com/noahbraunf/ysopy",
    packages=setuptools.find_packages(),
    classifiers=["Programming Language :: Python :: 3"],
    python_requires=">=3.6",
)

