DNA Curvature Analysis
======================

Dnacurve is a Python library and console script to calculate the global
3D structure of a B-DNA molecule from its nucleotide sequence according to the
dinucleotide wedge model. Local bending angles and macroscopic curvature
are calculated at each nucleotide.

For command line usage run ``python -m dnacurve --help``

:Author: `Christoph Gohlke <https://www.lfd.uci.edu/~gohlke/>`_

:License: BSD 3-Clause

:Version: 2021.6.29

Requirements
------------
* `CPython >= 3.7 <https://www.python.org>`_
* `Numpy 1.15.1 <https://www.numpy.org>`_
* `Matplotlib 3.2 <https://www.matplotlib.org>`_

Revisions
---------
2021.6.29
    Improve export to PDB.
2021.6.18
    Remove support for Python 3.6 (NEP 29).
    Fix dnacurve_web.py failure on WSL2.
2021.3.6
    Update copyright and formatting.
2020.1.1
    Remove support for Python 2.7 and 3.5.
    Update copyright.
2018.8.15
    Move modules into dnacurve package.
2018.5.29
    Add option to start web interface from console.
    Use matplotlib OOP interface.
2018.5.25
    Add functions to return PDB and CSV results as string.
2018.2.6
    Style and doctest fixes.
2014.6.16
    DNAse I Consensus model.
2013.11.21
    Overlapping chunks iterator.
2013.11.17
    Limit maximum sequence length to 510 nucleotides.
    Read simple Fasta sequence files.
    Save positive coordinates to PDB files.
    Fix sequence display for matplotlib 1.3.
2005.x.x
    Initial release.

Notes
-----
The algorithms, plots, and PDB format are not meant to be used with very
long sequences. By default sequences are truncated to 510 nucleotides,
which can be overridden by the user.

The generated PDB files can be visualized interactively using
`UCSF Chimera <https://www.cgl.ucsf.edu/chimera/>`_.

References
----------
1. Bending and curvature calculations in B-DNA.
   Goodsell DS, Dickerson RE. Nucleic Acids Res 22, 5497-503, 1994.
   See also http://mgl.scripps.edu/people/goodsell/research/bend/
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
>>> result.save_csv('_test.csv')
>>> result.save_pdb('_test.pdb')
>>> result.plot('_test.png', dpi=120)
