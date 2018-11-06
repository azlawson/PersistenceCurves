import setuptools

with open("README.txt", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PersistenceCurves",
    version="0.0.1",
    author="Austin Lawson",
    author_email="azlawson@uncg.edu",
    description="A small package created to aid in the calculation of Persistence Curves",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)