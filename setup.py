import setuptools

long_description = r'''
# celluloid
Effective Clustering for Single Cell Sequencing Cancer Data.
```
usage: celluloid [-h] {convert,cluster} ...

optional arguments:
  -h, --help         show this help message and exit

subcommands:
  valid subcommands

  {convert,cluster}
```
A detailed description of the module is available on our [github repo](https://github.com/AlgoLab/celluloid).
'''

setuptools.setup(
    name="celluloidclust",
    version="1.2.5",
    author="Murray Patterson",
    author_email="mpatterson@cs.gsu.edu",
    description="Effective Clustering for Single Cell Sequencing Cancer Data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AlgoLab/celluloid",
    packages=setuptools.find_packages(),
    install_requires=[
        'TatSu>=5.5.0',
        'kmodes>=0.10.2'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points={'console_scripts': [
        'celluloid=celluloid.celluloid:main']},
)