from setuptools import setup, find_packages
from pathlib import Path

def read_requirements(path = "requirements.txt"):
    p = Path(path)
    if not p.exists(): return []
    stripped_lines = (ln.strip() for ln in p.read_text().strip().splitlines())
    return [ln for ln in stripped_lines if ln and not ln.startswith("#")]

setup(
    name="volgrids",
    version="0.2.0",
    description="Framework for volumetric calculations, with emphasis in biological molecular systems.",
    keywords="grid mif smif volumetric molecular structural biology interaction field",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="DiegoBarMor",
    author_email="diegobarmor42@gmail.com",
    url="https://github.com/diegobarmor/volgrids",
    license="MIT",
    packages=find_packages(),
    package_data={"volgrids": ["_tables/*", "utils/*"]},
    install_requires=read_requirements(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)
