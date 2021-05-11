# from distutils.core import setup
from setuptools import setup, find_packages


# with open("README.md", "r") as fh:
#     long_description = fh.read()

setup(
    name="parablade",
    version="1.0",
    url='NA',
    author="NAnand-TUD&RoberAgro Roberto Agromayor",
    author_email="gingerxjiang@gmail.com",
    description="An open-source Python library for the parametrization of turbomachinery blades design",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
        ],
    packages=['parablade'], 
    zip_safe=False,   
    # py_modules=[
    #     "src"
    #     ]
    )