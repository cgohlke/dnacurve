DNA curvature analysis
======================

Dnacurve is a Python library, console script, and web application to calculate
the global 3D structure of a B-DNA molecule from its nucleotide sequence
according to the dinucleotide wedge model. Local bending angles and macroscopic
curvature are calculated at each nucleotide.

:Author: `Christoph Gohlke <https://www.cgohlke.com>`_
:License: BSD 3-Clause
:Version: 2023.8.30
:DOI: `10.5281/zenodo.7135499 <https://doi.org/10.5281/zenodo.7135499>`_

Quickstart
----------

Install the dnacurve package and all dependencies from the
`Python Package Index <https://pypi.org/project/dnacurve/>`_::

    python -m pip install -U dnacurve[all]

Print the console script usage::

    python -m dnacurve --help

Run the web application::

    python -m dnacurve --web

See `Examples`_ for using the programming interface.

Source code and support are available on
`GitHub <https://github.com/cgohlke/dnacurve>`_.

Requirements
------------

This revision was tested with the following requirements and dependencies
(other versions may work):

- `CPython <https://www.python.org>`_ 3.9.13, 3.10.11, 3.11.5, 3.12rc
- `Numpy <https://pypi.org/project/numpy/>`_ 1.25.2
- `Matplotlib <https://pypi.org/project/matplotlib/>`_ 3.7.2
- `Flask <https://pypi.org/project/Flask/>`_ 2.3.3 (optional)

Revisions
---------

2023.8.30

- Fix linting issues.
- Add py.typed marker.

2023.4.30

- Improve type hints.
- Drop support for Python 3.8 and numpy < 1.21 (NEP29).

2022.10.4

- Rename dnacurve_web.py to web.py (breaking).
- Deprecate save functions (use write functions).
- Add options to specify URL of web application and not opening web browser.
- Run web application using Flask if installed.
- Convert to Google style docstrings.
- Add type hints.
- Remove support for Python 3.7 and numpy < 1.19 (NEP29).

2021.6.29

- Improve export to PDB.

2021.6.18

- Remove support for Python 3.6 (NEP 29).
- Fix dnacurve_web.py failure on WSL2.

2021.3.6

- Update copyright and formatting.

2020.1.1

- Remove support for Python 2.7 and 3.5.
- Update copyright.

2018.8.15

- Move modules into dnacurve package.

2018.5.29

- Add option to start web interface from console.
- Use matplotlib OOP interface.

2018.5.25

- Add functions to return PDB and CSV results as string.

2018.2.6

- Style and doctest fixes.

2014.6.16

- DNAse I Consensus model.

2013.11.21

- Overlapping chunks iterator.

2013.11.17

- Limit maximum sequence length to 510 nucleotides.
- Read simple FASTA sequence files.
- Save positive coordinates to PDB files.
- Fix sequence display for matplotlib 1.3.

2005.x.x

- Initial release.

Notes
-----

The algorithms, plots, and PDB format are not meant to be used with very
long sequences. By default, sequences are truncated to 510 nucleotides,
which can be overridden by the user.

The generated PDB files can be visualized interactively using
`UCSF Chimera <https://www.cgl.ucsf.edu/chimera/>`_.

Dnacurve.py was derived from DNACG.PAS (c) 1993, and DNACURVE.CPP (c) 1995.

References
----------

1. Bending and curvature calculations in B-DNA.
   Goodsell DS, Dickerson RE. Nucleic Acids Res 22, 5497-503, 1994.
   See also http://mgl.scripps.edu/people/goodsell/research/bend/index.html.
2. Curved DNA without A-A: experimental estimation of all 16 DNA wedge angles.
   Bolshoy A et al. Proc Natl Acad Sci USA 88, 2312-6, 1991.
3. A comparison of six DNA bending models.
   Tan RK and Harvey SC. J Biomol Struct Dyn 5, 497-512, 1987.
4. Curved DNA: design, synthesis, and circularization.
   Ulanovsky L et al. Proc Natl Acad Sci USA 83, 862-6, 1986.
5. The ten helical twist angles of B-DNA.
   Kabsch W, Sander C, and Trifonov EN. Nucleic Acids Res 10, 1097-1104, 1982.
6. Rod models of DNA: sequence-dependent anisotropic elastic modelling of
   local bending phenomena.
   Munteanu MG et al. Trends Biochem Sci 23(9), 341-7, 1998.

Examples
--------

>>> from dnacurve import CurvedDNA
>>> result = CurvedDNA('ATGCAAATTG' * 5, 'trifonov', name='Example')
>>> result.curvature[:, 18:22]
array([[0.58062, 0.58163, 0.58278, 0.58378],
       [0.0803 , 0.11293, 0.07676, 0.03166],
       [0.57924, 0.5758 , 0.57368, 0.5735 ]])
>>> result.write_csv('_test.csv')
>>> result.write_pdb('_test.pdb')
>>> result.plot('_test.png', dpi=120)
