[![Build Status](https://travis-ci.org/NOAA-GFDL/MOM6.svg?branch=dev/master)](https://travis-ci.org/NOAA-GFDL/MOM6)
[![Read The Docs Status](https://readthedocs.org/projects/mom6/badge/?badge=latest)](http://mom6.readthedocs.io/)
[![codecov](https://codecov.io/gh/NOAA-GFDL/MOM6/branch/dev%2Fmaster/graph/badge.svg)](https://codecov.io/gh/NOAA-GFDL/MOM6)

# MOM6

This is the MOM6 source code used for the paper:

> Hailu Kong and Malte F. Jansen, 2021: Time-dependent response of the overturning circulation and pycnocline depth to Southern Ocean wind stress changes

that has been submitted to Journal of Physical Oceanography. Some modules are behind [NOAA-GFDL/MOM6](https://github.com/NOAA-GFDL/MOM6), where the present repository was forked from, but to ensure the reproducibility of simulations presented in the paper, which were done back in 2020 when MOM6 was still under frequent updates, one is recommended not to update the present repo with the most recent MOM6 version. 

# Where to find information

Start at the [MOM6-examples wiki](https://github.com/NOAA-GFDL/MOM6-examples/wiki) which has installation instructions.

For more specific installation instructions that are suitable to [Midway](https://rcc.uchicago.edu/docs/using-midway/index.html) super computing center at the University of Chicago, refer to [hlkong/MOM6-FMS-SIS2/wiki](https://github.com/hlkong/MOM6-FMS-SIS2/wiki). 

[Source code documentation](http://mom6.readthedocs.io/) is hosted on read the docs.

# What files are what

The top level directory structure groups source code and input files as follow:

| File/directory    | Purpose |
| --------------    | ------- |
| ```LICENSE.md```  | A copy of the Gnu lesser general public license, version 3. |
| ```README.md```   | This file with basic pointers to more information. |
| ```src/```        | Contains the source code for MOM6 that is always compiled. |
| ```config_src/``` | Contains optional source code depending on mode and configuration such as dynamic-memory versus static, ocean-only versus coupled. |
| ```pkg/```        | Contains third party (non-MOM6 or FMS) code that is compiled into MOM6. |
| ```docs/```       | Workspace for generated documentation. |

# Disclaimer

The United States Department of Commerce (DOC) GitHub project code is provided
on an "as is" basis and the user assumes responsibility for its use. DOC has
relinquished control of the information and no longer has responsibility to
protect the integrity, confidentiality, or availability of the information. Any
claims against the Department of Commerce stemming from the use of its GitHub
project will be governed by all applicable Federal law. Any reference to
specific commercial products, processes, or services by service mark,
trademark, manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce. The
Department of Commerce seal and logo, or the seal and logo of a DOC bureau,
shall not be used in any manner to imply endorsement of any commercial product
or activity by DOC or the United States Government.

This project code is made available through GitHub but is managed by NOAA-GFDL
at https://gitlab.gfdl.noaa.gov.
