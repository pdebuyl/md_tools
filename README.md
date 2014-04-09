md_tools: A few tools for Molecular Dynamics analysis
=====================================================

Copyright © 2013-2014 Pierre de Buyl

md_tools is a Python module providing a few computational routines to analyse
Molecular Dynamics trajectories.

md_tools is developped by Pierre de Buyl and is released under the modified BSD
license that can be found in the file LICENSE.

md_tools contains code from nMOLDYN, an interactive analysis program for
Molecular Dynamics simulations. It is especially designed for the computation
and decomposition of neutron scattering spectra, but also computes other
quantities.
nMOLDYN home page: http://dirac.cnrs-orleans.fr/nMOLDYN/
see licenses/LICENSE_nMOLDYN-3.0.10.txt file

md_tools documentation system is taken from lua-units by Peter Colberg
lua-units home page: http://colberg.org/lua-units/
see licenses/LICENSE_lua-units.txt file

Install
-------

md_tools requirements:
  - NumPy
  - Scientific version >= 2.8
  - Cython

Installation instructions

    python setup.py install --user

installs md_tools for the current user
