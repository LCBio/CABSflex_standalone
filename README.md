# CABS-flex

Welcome to the new **CABS-flex** in Python 3!

1. git clone [https://github.com/LCBio/cabsflex.git](https://github.com/LCBio/cabsflex/tree/main)
2. conda env create -f environment.yml
3. conda activate flex_prod
4. pip install .
5. Download the [NetSurfP-3.0_standalone](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netsurfp&version=3.0&packageversion=3.0&platform=Linux).
6. Unzip the downloaded file and copy the `NetSurfP-3.0_standalone` folder to the `cabsflex` folder.
7. cd ./NetSurfP-3.0_standalone
8. python setup.py install