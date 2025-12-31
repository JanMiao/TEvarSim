from setuptools import setup, find_packages

setup(
    name="TEvarSim",
    version="1.0.3",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy",
        "biopython",
        "pysam"
    ],
    entry_points={
        "console_scripts": [
            "tevarsim=TEvarSim.__main__:main"
        ]
    },
    description="TEvarSim is a versatile genome simulation tool for generating polymorphic transposable element (TE) variants.",
    long_description="TEvarSim is a toolkit to simulate and analyze polymorphic transposable elements across genomes, supporting TE pool building, real TE processing, pangenome TE extraction, pTE simulation, and VCF comparison.",
    long_description_content_type="text/markdown",
    author="JIAN MIAO",
    author_email="miaojian6363@gmail.com",
    url="https://github.com/JanMiao/TEvarSim",  
    classifiers=[               
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
