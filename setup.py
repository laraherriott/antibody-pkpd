#
# setuptools script
#
from setuptools import setup, find_packages

setup(
    # Module name (lowercase)
    name='lecanemab_model',

    description='Implemetation of a lecanemba PK/PD model.',

    license='MIT license',

    # author='',

    # author_email='',

    maintainer='Lara Herriott',

    maintainer_email='lara.herriott@dtc.ox.ac.uk',

    # Packages to include
    packages=find_packages(include=('lecanemab', 'lecanemab.*')),

    # List of dependencies
    install_requires=[
        # Dependencies go here!
        'numpy',
        'matplotlib',
        'scipy',
    ]
)