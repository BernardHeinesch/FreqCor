from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="FreqCor",
    version="1.3",
    author="Ariane Faurès",
    author_email="bernard.heinesch@uliege.be",
    description="A Python toolbox to compute and apply spectral/cospectral correction factors to eddy covariance fluxes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BernardHeinesch/FreqCor",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
)
