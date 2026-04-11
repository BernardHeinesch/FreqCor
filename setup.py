from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="FreqCor",
    version="1.4",
    author="Ariane Faurès",
    author_email="bernard.heinesch@uliege.be",
    description="A Python toolbox to compute and apply spectral/cospectral correction factors to eddy covariance fluxes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BernardHeinesch/FreqCor",
    py_modules=[
        "FREQCOR_cof",
        "FREQCOR_Compute",
        "FREQCOR_Flux",
        "FREQCOR_functions",
        "FREQCOR_LUT_CF",
        "FREQCOR_LUT_cof",
        "FREQCOR_Main",
        "FREQCOR_plot",
        "FREQCOR_Read_EP",
        "FREQCOR_Read_GEddySoft",
        "FREQCOR_Read_TOF",
        "FREQCOR_Sel_CF",
        "FREQCOR_Sel_cof",
        "FREQCOR_Sel_general",
        "FREQCOR_Sel_stunst",
        "FREQCOR_Sensor_Separation",
        "FREQCOR_Start",
        "FREQCOR_validate",
        "FREQCOR_VM_flag",
        "FREQCOR_write_outputs",
        "FREQCOR_Ref_cospectrum_for_plotting",
    ],
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

