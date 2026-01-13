from setuptools import setup, find_packages

setup(
    name="volgrids",
    version="0.1.0",
    description="Framework for volumetric calculations, with emphasis in biological molecular systems.",
    keywords="grid",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="DiegoBarMor",
    author_email="diegobarmor42@gmail.com",
    url="https://github.com/diegobarmor/volgrids",
    license="MIT",
    packages=find_packages(),
    package_data={
        "volgrids": ["_data/*"],
    },
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)
