import setuptools

with open("README.txt", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="persistencecurves",
    version="0.0.7",
    author="Austin Lawson",
    author_email="azlawson@uncg.edu",
    description="A small package created to aid in the calculation of Persistence Curves",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    url="https://github.com/azlawson/PersistenceCurves",
    install_requires = ["numpy", "scipy","matplotlib"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)