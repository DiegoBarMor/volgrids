from setuptools import setup, find_packages
from pathlib import Path

__version__: str
exec(Path("volgrids/_version.py").read_text())

setup(
    name="volgrids",
    version=__version__,
    description="Framework for volumetric calculations, with emphasis in biological molecular systems.",
    keywords="grid mif smif volumetric molecular structural biology interaction field",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="DiegoBarMor",
    author_email="diegobarmor42@gmail.com",
    url="https://github.com/diegobarmor/volgrids",
    license="MIT",
    packages=find_packages(),
    package_data={"volgrids": ["_tables/*", "apbs/*"]},
    install_requires=["MDAnalysis==2.10.0", "h5py==3.15.1"],
    entry_points={
        "console_scripts": [
            "volgrids=volgrids.__main__:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)
