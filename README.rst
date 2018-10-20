DNA Curvature Analysis
======================

Dnacurve is a Python library and console script to calculate the global
3D structure of a B-DNA molecule from its nucleotide sequence according to the
dinucleotide wedge model. Local bending angles and macroscopic curvature
are calculated at each nucleotide.

For command line usage run ``python -m dnacurve --help``

:Author: `Christoph Gohlke <https://www.lfd.uci.edu/~gohlke/>`_

:Version: 2018.10.18

Requirements
------------
* `CPython 2.7 or 3.5+ <https://www.python.org>`_
* `Numpy 1.14 <https://www.numpy.org>`_
* `Matplotlib 2.2 <https://www.matplotlib.org>`_

Revisions
---------
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
The API is not stable yet and is expected to change between revisions.

The algorithms, plots, and PDB format are not meant to be used with very
long sequences. By default sequences are truncated to 510 nucleotides,
which can be overridden by the user.

The generated PDB files can be visualized interactively using
`UCSF Chimera <http://www.cgl.ucsf.edu/chimera/>`_.

References
----------
(1) Bending and curvature calculations in B-DNA.
    Goodsell DS, Dickerson RE. Nucleic Acids Res 22, 5497-503, 1994.
    See also http://mgl.scripps.edu/people/goodsell/research/bend/
(2) Curved DNA without A-A: experimental estimation of all 16 DNA wedge angles.
    Bolshoy A et al. Proc Natl Acad Sci USA 88, 2312-6, 1991.
(3) A comparison of six DNA bending models.
    Tan RK and Harvey SC. J Biomol Struct Dyn 5, 497-512, 1987.
(4) Curved DNA: design, synthesis, and circularization.
    Ulanovsky L et al. Proc Natl Acad Sci USA 83, 862-6, 1986.
(5) The ten helical twist angles of B-DNA.
    Kabsch W, Sander C, and Trifonov EN. Nucleic Acids Res 10, 1097-1104, 1982.
(6) Rod models of DNA: sequence-dependent anisotropic elastic modelling of
    local bending phenomena.
    Munteanu MG et al. Trends Biochem Sci 23(9), 341-7, 1998.

Examples
--------
>>> from dnacurve import CurvedDNA
>>> result = CurvedDNA('ATGCAAATTG'*5, 'trifonov', name='Example')
>>> result.curvature[:, 18:22]
array([[ 0.58061616,  0.58163338,  0.58277938,  0.583783  ],
       [ 0.08029914,  0.11292516,  0.07675816,  0.03166286],
       [ 0.57923902,  0.57580064,  0.57367815,  0.57349872]])
>>> result.save_csv('_test.csv')
>>> result.save_pdb('_test.pdb')
>>> result.plot('_test.png', dpi=160)
