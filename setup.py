from setuptools import setup
from CABS import __version__

setup(
    name='CABS',
    version=__version__,
    packages=['CABS'],
    url='https://github.com/LCBio/cabsflex',
    license='MIT',
    author='Laboratory of Computational Biology',
    author_email='k.wroblewski7@uw.edu.pl',
    description='CABS in python3',
    entry_points={
        'console_scripts': [
            'CABSdock = CABS.__main__:run_dock',
            'CABSflex = CABS.__main__:run_flex'
        ]
    },
    package_data={'CABS': ['data/*.dat', 'config.json']}
)
