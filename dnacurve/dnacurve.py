# dnacurve.py

# Copyright (c) 1993-2024, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""DNA curvature analysis.

Dnacurve is a Python library, console script, and web application to calculate
the global 3D structure of a B-DNA molecule from its nucleotide sequence
according to the dinucleotide wedge model. Local bending angles and macroscopic
curvature are calculated at each nucleotide.

:Author: `Christoph Gohlke <https://www.cgohlke.com>`_
:License: BSD 3-Clause
:Version: 2024.5.24
:DOI: `10.5281/zenodo.7135499 <https://doi.org/10.5281/zenodo.7135499>`_

Quickstart
----------

Install the dnacurve package and all dependencies from the
`Python Package Index <https://pypi.org/project/dnacurve/>`_::

    python -m pip install -U "dnacurve[all]"

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

- `CPython <https://www.python.org>`_ 3.9.13, 3.10.11, 3.11.9, 3.12.3
- `Numpy <https://pypi.org/project/numpy/>`_ 1.26.4
- `Matplotlib <https://pypi.org/project/matplotlib/>`_ 3.8.4
- `Flask <https://pypi.org/project/Flask/>`_ 3.0.3 (optional)

Revisions
---------

2024.5.24

- Fix docstring examples not correctly rendered on GitHub.

2024.5.10

- Fix mypy errors.

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
>>> cdna = CurvedDNA('ATGCAAATTG' * 5, 'trifonov', name='Example')
>>> cdna.curvature[:, 18:22]
array([[0.58062, 0.58163, 0.58278, 0.58378],
       [0.0803 , 0.11293, 0.07676, 0.03166],
       [0.57924, 0.5758 , 0.57368, 0.5735 ]])
>>> cdna.write_csv('_test.csv')
>>> cdna.write_pdb('_test.pdb')
>>> cdna.plot('_test.png', dpi=120)

"""

from __future__ import annotations

__version__ = '2024.5.24'

__all__ = [
    'CurvedDNA',
    'Model',
    'Sequence',
    'MAXLEN',
    'MODELS',
    'chunks',
    'complementary',
    'dinuc_window',
    'dinucleotide_matrix',
    'main',
    'oligonucleotides',
    'overlapping_chunks',
    'superimpose_matrix',
    'unique_oligos',
]

import datetime
import math
import os
import re
import sys
import warnings
from typing import TYPE_CHECKING, overload

import numpy

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from typing import Any, BinaryIO

    from numpy.typing import NDArray

MAXLEN: int = 510
"""Maximum length of sequences to analyze."""

MODELS: list[str] = []  # updated later
"""Predefined models."""


class CurvedDNA:
    """Calculate, plot or write helix coordinates, local bending and curvature.

    Parameters:
        sequence:
            Sequence instance, file name, or nucleotide sequence.
            See :py:class:`Sequence` parameters.
        model:
            Model instance, file name, class, dict, or name of
            predefined model. See :py:class:`Model` parameters.
        name:
            Name of sequence.
        curvature_window:
            Window size for calculating curvature.
        bend_window:
            Window size for calculating local bend angles.
        curve_window:
            Window size for calculating curvature angles.
        maxlen:
            Maximum length of sequence.

    Notes:
        Atomic coordinates are centered at origin and oriented such that:
        (1) helix-axis endpoints lie on x-axis and
        (2) maximum deviation of DNA- from x-axis is along the z-axis.
        Coordinates in PDB files are shifted to the positive domain.

        The **curvature** at nucleotide N is one over the radius of a
        circle passing through helix axis coordinates N-window, N, and
        N+window, which are separated by one respectively two helix turns.
        The three points define a triangle.  The radius is the product of
        the length of the triangle sides divided by four times the area of
        the triangle.  A window size of 10 is optimal for B-DNA.

        The **local bend angle** at nucleotide N is the angle between the
        normal vectors of base pairs N-window and N+window.  The window size
        should be one or two.

        The **curvature angle** at nucleotide N is the angle between the
        smoothed normal vectors of basepair N-window and N+window.
        The window size should be in the order of 15.

        The curvature and bend values are normalized relative to the
        DNA curvature in a nucleosome (0.0234).

    """

    P_COORDINATES = (
        8.91,
        -5.2,
        2.08,
    )
    """Cylindrical coordinates of 5' phosphate:

    0. distance from axis
    1. angle to roll axis
    2. distance from basepair plane
    """

    coordinates: NDArray[Any]
    """Homogeneous coordinates at each nucleotide of:

    0. helix axis
    1. phosphate of 5'-3' strand
    2. phosphate of antiparallel strand
    3. basepair normal vector
    4. smoothed basepair normal vector
    """

    curvature: NDArray[Any]
    """Values at each nucleotide relative to curvature in nucleosome:

    0. curvature
    1. local bend angle
    2. curvature angle
    """

    scales: NDArray[Any]
    """Scaling factors used to normalize curvature array."""

    windows: tuple[int, int, int]
    """Window sizes for calculating curvature, local bend angle, and
    curvature angle.
    """

    date: datetime.datetime
    """Date and time of analysis."""

    _limits: tuple[float, float, float]

    def __init__(
        self,
        sequence: Sequence | os.PathLike[Any] | str,
        /,
        model: Model | os.PathLike[Any] | str = 'trifonov',
        name: str = 'Untitled',
        curvature_window: int = 10,
        bend_window: int = 2,
        curve_window: int = 15,
        maxlen: int = MAXLEN,
    ) -> None:
        if isinstance(model, Model):
            self.model = model
        else:
            self.model = Model(model)
        if isinstance(sequence, Sequence):
            self.sequence = sequence
        else:
            self.sequence = Sequence(sequence, name, maxlen=maxlen)
        size = len(self.sequence)
        if size < self.model.order:
            raise ValueError(
                f'sequence must be >{self.model.order} nucleotides long'
            )
        if size > maxlen:
            warnings.warn(f'sequence is longer than {maxlen} nucleotides')

        assert 0 < curvature_window < 21
        assert 0 < bend_window < 4
        assert 9 < curve_window < 21
        self.windows = (curvature_window, bend_window, curve_window)
        self._limits = (10.0, 10.0, 10.0)

        self.date = datetime.datetime.now()
        self.coordinates = numpy.zeros((5, size, 4), dtype=numpy.float64)
        self.curvature = numpy.zeros((3, size), dtype=numpy.float64)
        self.scales = numpy.ones((3, 1), dtype=numpy.float64)

        self._coordinates()
        self._reorient()
        self._center()
        self._curvature()

    def _coordinates(self) -> None:
        """Calculate coordinates and normal vectors from sequence and model."""
        s, a, z = self.P_COORDINATES
        x = s * math.cos(math.radians(a))
        y = s * math.sin(math.radians(a))
        xyz = self.coordinates
        xyz[0:3, :, 3] = 1.0  # homogeneous coordinates
        xyz[1, :, 0:3] = x, y, z  # 5' phosphate
        xyz[2, :, 0:3] = -x, y, -z  # phosphate of antiparallel strand
        xyz[3, :, 2] = 1.0  # basepair normal vectors

        matrices = self.model.matrices
        dot = numpy.dot
        for i, seq in enumerate(dinuc_window(self.sequence, self.model.order)):
            xyz[:4, : i + 1, :] = dot(xyz[:4, : i + 1, :], matrices[seq])

        # average direction vector of one helix turn,
        # calculated by smoothing the basepair normals
        if len(self.sequence) > 10:
            kernel = numpy.array([0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5])
            kernel /= kernel.sum()
            for i in 0, 1, 2:
                xyz[4, :, i] = numpy.convolve(kernel, xyz[3, :, i], 'same')
            for i in range(5, len(self.sequence) - 5):
                xyz[4, i, :] /= norm(xyz[4, i, :])

    def _reorient(self) -> None:
        """Reorient coordinates."""
        xyz: NDArray[Any] = self.coordinates[0, :, 0:3]  # helix axis
        xyz = xyz - xyz[-1]
        # assert start point is at origin
        assert numpy.allclose(xyz[-1], (0, 0, 0))
        # normalized end to end vector
        e: NDArray[Any] = +xyz[0]
        e_len = norm(e)
        e /= e_len
        # point i of maximum distance to end to end line
        t: NDArray[Any] = numpy.cross(e, xyz)
        t = numpy.sum(t * t, axis=1)
        i = numpy.argmax(t)
        x = math.sqrt(t[i])
        # distance of endpoint to xyz[i]
        w = norm(xyz[i])
        # distance of endpoint to point on end to end line nearest to xyz[i]
        u = math.sqrt(w * w - x * x)
        # find transformation matrix
        v0: NDArray[Any] = xyz[[0, i, -1]]  # type: ignore
        v1: NDArray[Any] = numpy.array(
            ((0, 0, 0), (e_len - u, 0, x), (e_len, 0, 0))
        )
        M = superimpose_matrix(v0, v1)
        self.coordinates: NDArray[Any] = numpy.dot(self.coordinates, M.T)

    def _center(self) -> None:
        """Center atomic coordinates at origin."""
        xyz = self.coordinates[0:3, :, 0:3]  # helix axis and P atoms
        low = numpy.min(numpy.min(xyz, axis=1), axis=0)
        upp = numpy.max(numpy.max(xyz, axis=1), axis=0)
        lim = (upp - low) / 2.0
        self.coordinates[0:3, :, 0:3] -= low + lim
        self._limits = float(lim[0]), float(lim[1]), float(lim[2])

    def _curvature(self) -> None:
        """Calculate normalized curvature and bend angles."""
        dot = numpy.dot
        cross = numpy.cross
        arccos = numpy.arccos
        size = len(self.sequence)

        # curvature from radius
        window = self.windows[0]
        if size >= 2 * window:
            result = self.curvature[0, :]
            xyz = self.coordinates[0, :, 0:3]  # helix axis
            for i in range(window, size - window):
                a, b, c = xyz[[i - window, i, i + window]]
                ab = b - a
                bc = c - b
                ac = c - a
                lab = norm(ab)
                lbc = norm(bc)
                lac = norm(ac)
                area = norm(cross(ab, ac))
                result[i] = (2.0 * area) / (lab * lbc * lac)

        # local bend angles from basepair normals
        window = self.windows[1]
        if size >= 2 * window:
            normals = self.coordinates[3, :, 0:3]
            result = self.curvature[1, :]
            for i in range(window, size - window):
                if not numpy.allclose(
                    normals[i - window], normals[i + window]
                ):
                    result[i] = arccos(
                        dot(normals[i - window], normals[i + window])
                    )

        # curvature angles from smoothed basepair normals
        window = self.windows[2]
        if size >= 2 * window:
            normals = self.coordinates[4, :, 0:3]
            result = self.curvature[2, :]
            for i in range(window, size - window):
                if not numpy.allclose(
                    normals[i - window], normals[i + window]
                ):
                    result[i] = arccos(
                        dot(normals[i - window], normals[i + window])
                    )

        self.curvature = numpy.nan_to_num(self.curvature)

        # normalize relative to curvature in nucleosome
        self.scales[0] = 0.0234
        self.scales[1:] = 0.0234 * 2 * self.windows[2] * self.model.rise
        self.curvature /= self.scales

    def write_csv(self, path: os.PathLike[Any] | str, /) -> None:
        """Write coordinates and curvature values to CSV file.

        Parameters:
            path: Name of CSV file to write.

        """
        with open(path, 'w', encoding='latin-1', newline='\r\n') as fh:
            fh.write(self.csv())

    def save_csv(self, path: os.PathLike[Any] | str) -> None:
        # deprecated
        self.write_csv(path)

    def write_pdb(self, path: os.PathLike[Any] | str, /) -> None:
        """Write atomic coordinates to PDB file.

        Parameters:
            path: Name of PDB file to write.

        """
        with open(path, 'w', encoding='latin-1', newline='\n') as fh:
            fh.write(self.pdb())

    def save_pdb(self, path: os.PathLike[Any] | str) -> None:
        # deprecated
        self.write_pdb(path)

    def csv(self) -> str:
        """Return coordinates and curvature values in CSV format."""
        seq = self.sequence
        cur = self.curvature
        xyz = self.coordinates
        csv = [
            'Sequence,Index,Curvature / {:.4f} [{}],'
            'Bend Angle / {:.4f}  [{}],Curvature Angle / {:.4f}  [{}],'
            'Helix Axis [x],Helix Axis [y],Helix Axis [z],'
            'Phosphate 1 [x],Phosphate 1 [y],Phosphate 1 [z],'
            'Phosphate 2 [x],Phosphate 2 [y],Phosphate 2 [z],'
            'Basepair Normal [x],Basepair Normal [y],Basepair Normal [z],'
            'Smoothed Normal [x],Smoothed Normal [y],Smoothed Normal [z]'
            '\n'.format(
                self.scales[0, 0],
                self.windows[0],
                self.scales[1, 0],
                self.windows[1],
                self.scales[2, 0],
                self.windows[2],
            )
        ]
        for i in range(len(self.sequence)):
            csv.append(f'{seq[i]},{i + 1}')
            csv.append(',{:.5f},{:.5f},{:.5f}'.format(*cur[:, i]))
            for j in range(5):
                csv.append(',{:.5f},{:.5f},{:.5f}'.format(*xyz[j, i, 0:3]))
            csv.append('\n')
        return ''.join(csv)

    def pdb(self) -> str:
        """Return atomic coordinates in PDB format."""
        pdb: list[str] = []

        def pdb_append(line: str) -> None:
            pdb.append(f'{line:<80s}\n')

        date1 = datetime.datetime.strftime(self.date, '%d-%b-%y').upper()
        date2 = datetime.datetime.strftime(self.date, '%d %b %Y, %H:%M:%S')
        pdb_append(
            f'HEADER    CURVED DNA                              {date1}   0XYZ'
        )
        pdb_append('TITLE     RESULT OF DNA CURVATURE ANALYSIS')
        pdb_append('COMPND    MOL_ID: 1;')
        pdb_append('COMPND   2 CHAIN: A, B;')
        pdb_append('EXPDTA    THEORETICAL MODEL')
        pdb_append(f'REMARK   1 GENERATED BY DNACURVE.PY ON {date2.upper()}')
        pdb_append(f'REMARK   2 PARAMETERS: {self.model.name.upper()[:48]}')
        pdb_append(f'REMARK   3 SEQUENCE: {self.sequence.name.upper()[:48]}')
        seq = self.sequence
        seqi = complementary(seq)
        seqlen = len(seq)
        xyz = self.coordinates[0:3, :, 0:3]
        xyz = xyz - xyz.min(axis=1).min(axis=0)

        for i, chunk in enumerate(chunks(seq.string, 13)):
            residues = ' '.join(f' D{s}' for s in chunk)
            pdb_append(f'SEQRES {i + 1:>3} A {seqlen:>4}  {residues}')
        for i, chunk in enumerate(chunks(seqi, 13)):
            residues = ' '.join(f' D{s}' for s in chunk)
            pdb_append(f'SEQRES {i + 1:>3} B {seqlen:>4}  {residues}')

        serial = 1
        for i in range(seqlen):
            x, y, z = xyz[1, i, 0:3]
            pdb_append(
                f'ATOM  {serial:5}  P    D{seq[i]:<1s} A{serial:4}'
                f'    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           P'
            )
            serial += 1
        pdb_append(f'TER   {serial:5}       D{seq[i]:<1s} A{serial - 1:4}')

        serial += 1
        for i in range(seqlen):
            x, y, z = xyz[2, seqlen - i - 1, 0:3]
            pdb_append(
                f'ATOM  {serial:5}  P    D{seqi[i]:<1s} B{serial:4}'
                f'    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           P'
            )
            serial += 1
        pdb_append(f'TER   {serial:5}       D{seq[i]:<1s} B{serial - 1:4}')
        pdb_append('END')
        return ''.join(pdb)

    def plot(
        self,
        arg: str | os.PathLike[Any] | BinaryIO | bool = True,
        /,
        dpi: int = 96,
        figsize: tuple[float, float] = (6.0, 7.5),
        imageformat: str | None = None,
    ) -> Any:
        """Plot results using matplotlib.

        Parameters:
            arg:
                Specifies how to plot:

                - **False**: do not plot.
                - **True**: plot interactively.
                - **File name**: write the figure to it.
                - **Open file**: write the figure to it.
            dpi:
                Resolution of plot in dots per inch.
            figsize:
                Matplotlib figure size.
            imageformat:
                Image file format of the figure, such as, 'png', 'pdf',
                or 'svg'.

        """
        if not arg:
            return None

        if arg is True:
            # GUI plot
            from matplotlib import pyplot

            fig = pyplot.figure(dpi=dpi, figsize=figsize)
            fcm = fig.canvas.manager
            if fcm is not None:
                fcm.set_window_title('DNA Curvature Analysis')
        else:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            from matplotlib.figure import Figure

            fig = Figure(dpi=dpi, figsize=figsize)
            FigureCanvasAgg(fig)

        # sequence
        ax = fig.add_subplot(321)
        ax.set_title(
            f'Sequence: {self.sequence.name}\nModel: {self.model.name}',
            x=0,
            horizontalalignment='left',
            size=11,
        )
        ax.text(
            0.0,
            0.95,
            self.sequence.format(line=3),
            family='monospace',
            size=7,
            verticalalignment='top',
            horizontalalignment='left',
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        ax.axis('image')

        # projections
        xyz = self.coordinates[..., 0:3]
        limit = max(self._limits) * 1.1

        def plot_projection(
            plotnum: int, axes: str, label: bool = True, /
        ) -> None:
            ax = fig.add_subplot(plotnum)
            ax.set_title(f'{axes[0]}-{axes[1]}', size=11)
            ax0, ax1 = (ord(a) - 88 for a in axes)
            # TODO: sort lines by depth and draw back to front
            ax.plot(
                xyz[1, :, ax0],
                xyz[1, :, ax1],
                'r-',
                xyz[2, :, ax0],
                xyz[2, :, ax1],
                'b-',
                xyz[0, :, ax0],
                xyz[0, :, ax1],
                'k-',
                lw=0.8,
            )
            if label:
                ax.text(
                    xyz[1, 0, ax0],
                    xyz[1, 0, ax1],
                    "5'",
                    color='darkred',
                    clip_on=True,
                    size=10,
                )
                ax.text(
                    xyz[1, -1, ax0],
                    xyz[1, -1, ax1],
                    "3'",
                    color='darkred',
                    clip_on=True,
                    size=10,
                )
            ax.axis('image')
            ax.axis((-limit, limit, -limit, limit))
            ax.axis('off')

        plot_projection(322, 'XZ')
        plot_projection(323, 'ZY', False)
        plot_projection(324, 'XY')

        # plot
        ax = fig.add_subplot(325)
        ax.set_position((0.1, 0.075, 0.83, 0.21))

        def plot_curvature(
            index: int, label: str, style: str, lw: float, /
        ) -> None:
            w = self.windows[index]
            ax.plot(
                numpy.arange(w, len(self.sequence) - w),
                self.curvature[index, w:-w],
                style,
                lw=lw,
                label=label + f' / {self.scales[index, 0]:.4f} [{w}]',
            )

        plot_curvature(0, 'Curvature', 'r-', 0.8)
        plot_curvature(2, 'Curvature Angle', 'b-', 0.8)
        plot_curvature(1, 'Bend Angle', 'k-', 0.25)

        lg = ax.legend(loc=2, borderaxespad=0.0, prop={'size': 8})
        lg.get_frame().set_linewidth(0)
        lg.get_frame().set_fill(False)

        ax.get_xaxis().tick_bottom()
        ax.axhline(0.5, ls=':', color='0.5')
        ax.set_yticks((0, 0.5, 1))
        ax.tick_params(labelsize=8)
        ax.axis((0, len(self.sequence), 0, 1.0))
        ax.set_xlabel('Sequence Index', size=11)

        if arg is True:
            pyplot.show()
        elif arg:
            fig.savefig(arg, dpi=dpi, format=imageformat)
        else:
            return fig
        return None

    @property
    def name(self) -> str:
        """Name of sequence."""
        return self.sequence.name

    def __repr__(self) -> str:
        return (
            f'<{self.__class__.__name__} '
            f'{self.sequence.name!r} {self.model.name!r} '
            f'length={len(self.sequence)!r}>'
        )

    def __str__(self) -> str:
        return f'{self.sequence}\n\n{self.model}\n'


class Model:
    """N-mer DNA-bending model.

    Transformation parameters and matrices for all oligonucleotides of
    certain length.

    Parameters:
        model:
            Specifies type of model:

            - **None**: :py:attr:`Model.STRAIGHT`.
            - **Model**: :py:class:`Model` instance.
            - **dict**: Specifies :py:attr:`Model.name`, :py:attr:`Model.oligo`
              :py:attr:`Model.twist`, :py:attr:`Model.roll`,
              :py:attr:`Model.tilt`, and :py:attr:`Model.rise`.
            - **str**: Name of predefined model, 'straight', 'aawedge',
              'trifonov', 'desantis', 'calladine', or 'reversed'.
            - **Path**: File containing model parameters.
        **kwargs:
            Additional arguments to modify model:

            - :py:attr:`Model.name`
            - :py:attr:`Model.oligo`
            - :py:attr:`Model.twist`
            - :py:attr:`Model.roll`
            - :py:attr:`Model.tilt`
            - :py:attr:`Model.rise`

    Examples:
        >>> m = Model('AAWedge')
        >>> m = Model('Nucleosome')
        >>> m = Model(**Model.STRAIGHT)
        >>> m = Model(Model.CALLADINE, name='My Model', rise=4.0)
        >>> m.name == 'My Model' and m.rise == 4.0
        True
        >>> m = Model(
        ...     name='Test',
        ...     rise=3.38,
        ...     oligo='AA AC AG AT CA GG CG GA GC TA'.split(),
        ...     twist=(34.29,) * 10,
        ...     roll=(0.0,) * 10,
        ...     tilt=(0.0,) * 10,
        ... )
        >>> m.write('_test.dat')
        >>> m.twist == Model('_test.dat').twist
        True

    """

    # fmt: off
    # noqa: E201,E222,E241,E251
    STRAIGHT: dict[str, Any] = dict(
        name = 'Straight',
        oligo = 'AA AC AG AT CA GG CG GA GC TA',
        twist = (360.0 / 10.5, ) * 10,
        roll =  (0.0, ) * 10,
        tilt =  (0.0, ) * 10,
        rise = 3.38,
    )

    AAWEDGE: dict[str, Any] = dict(
        name = 'AA Wedge',
        oligo = 'AA     AC    AG    AT    CA    GG     CG    GA    GC    TA',
        twist = (35.62, 34.4, 27.7, 31.5, 34.5, 33.67, 29.8, 36.9, 40.0, 36.0),
        roll =  (-8.40,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0),
        tilt =  ( 2.40,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0),
        rise = 3.38,
    )

    CALLADINE: dict[str, Any] = dict(
        # Nucleic Acids Res, 1994, 22(24), p 5498, Table 1, Model b.
        name = 'Calladine & Drew',
        oligo = 'AA    AC    AG    AT    CA    GG    CG    GA    GC    TA',
        twist = (35.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0),
        roll =  ( 0.0,  3.3,  3.3,  3.3,  3.3,  3.3,  3.3,  3.3,  3.3,  6.6),
        tilt =  ( 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0),
        rise = 3.38,
    )

    TRIFONOV: dict[str, Any] = dict(
        # Nucleic Acids Res, 1994, 22(24), p 5498, Table 1, Model c.
        name = 'Bolshoi & Trifonov',
        oligo = 'AA     AC    AG    AT    CA    GG     CG    GA    GC    TA',
        twist = (35.62, 34.4, 27.7, 31.5, 34.5, 33.67, 29.8, 36.9, 40.0, 36.0),
        roll =  (-6.50, -0.9,  8.4,  2.6,  1.6,   1.2,  6.7, -2.7, -5.0,  0.9),
        tilt =  ( 3.20, -0.7, -0.3,  0.0,  3.1,  -1.8,  0.0, -4.6,  0.0,  0.0),
        rise = 3.38,
    )

    DESANTIS: dict[str, Any] = dict(
        # Nucleic Acids Res, 1994, 22(24), p 5498, Table 1, Model d.
        name = 'Cacchione & De Santis',
        oligo = 'AA    AC    AG    AT    CA    GG    CG    GA    GC    TA',
        twist = (35.9, 34.6, 35.6, 35.0, 34.5, 33.0, 33.7, 35.8, 33.3, 34.6),
        roll =  (-5.4, -2.4,  1.0, -7.3,  6.7,  1.3,  4.6,  2.0, -3.7,  8.0),
        tilt =  (-0.5, -2.7, -1.6,  0.0,  0.4, -0.6,  0.0, -1.7,  0.0,  0.0),
        rise = 3.38,
    )

    REVERSED: dict[str, Any] = dict(
        # Nucleic Acids Res, 1994, 22(24), p 5498, Table 1, Model e.
        name = 'Reversed Calladine & Drew',
        oligo = 'AA    AC    AG    AT    CA    GG    CG    GA    GC    TA',
        twist = (35.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0),
        roll =  ( 3.3,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -3.3),
        tilt =  ( 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0),
        rise = 3.38,
    )

    NUCLEOSOME: dict[str, Any] = dict(
        # Nucleic Acids Res, 1994, 22(24), p 5498, Table 1, Model a.
        name = 'Nucleosome Positioning',
        oligo = """
                AAA  ATA   AGA   ACA  TAA  TTA  TGA  TCA
                GAA  GTA   GGA   GCA  CAA  CTA  CGA  CCA
                AAT  ATT   AGT   ACT  TAT  TTT  TGT  TCT
                GAT  GTT   GGT   GCT  CAT  CTT  CGT  CCT
                AAG  ATG   AGG   ACG  TAG  TTG  TGG  TCG
                GAG  GTG   GGG   GCG  CAG  CTG  CGG  CCG
                AAC  ATC   AGC   ACC  TAC  TTC  TGC  TCC
                GAC  GTC   GGC   GCC  CAC  CTC  CGC  CCC""",
        roll = (0.0, 2.8,  3.3,  5.2, 2.0, 2.0, 5.4, 5.4,
                3.0, 3.7,  3.8,  5.4, 3.3, 2.2, 8.3, 5.4,
                0.7, 0.7,  5.8,  5.8, 2.8, 0.0, 5.2, 3.3,
                5.3, 3.7,  5.4,  7.5, 6.7, 5.2, 5.4, 5.4,
                5.2, 6.7,  5.4,  5.4, 2.2, 3.3, 5.4, 8.3,
                5.4, 6.5,  6.0,  7.5, 4.2, 4.2, 4.7, 4.7,
                3.7, 5.3,  7.5,  5.4, 3.7, 3.0, 6.0, 3.8,
                5.4, 5.4, 10.0, 10.0, 6.5, 5.4, 7.5, 6.0),
        twist = (34.3, ) * 64,
        tilt =  (0.0, ) * 64,
        rise = 3.38,
    )

    TRINUCLEOTIDE: dict[str, Any] = dict(
        # Trends Biochem Sci, 1998, 23(9), p 341, Table 1, consensus scale.
        name = 'DNAse I Consensus',
        oligo = """
                AAA   ATA   AGA   ACA   TAA   TTA   TGA   TCA
                GAA   GTA   GGA   GCA   CAA   CTA   CGA   CCA
                AAT   ATT   AGT   ACT   TAT   TTT   TGT   TCT
                GAT   GTT   GGT   GCT   CAT   CTT   CGT   CCT
                AAG   ATG   AGG   ACG   TAG   TTG   TGG   TCG
                GAG   GTG   GGG   GCG   CAG   CTG   CGG   CCG
                AAC   ATC   AGC   ACC   TAC   TTC   TGC   TCC
                GAC   GTC   GGC   GCC   CAC   CTC   CGC   CCC""",
        roll = (0.05, 6.25, 4.90, 5.50, 4.65, 4.65, 7.70, 7.70,
                4.05, 5.05, 5.00, 6.75, 4.75, 5.00, 7.05, 3.05,
                0.35, 0.35, 3.90, 3.90, 6.25, 0.05, 5.50, 4.90,
                4.45, 2.65, 5.30, 6.90, 7.70, 4.70, 5.30, 5.05,
                4.70, 7.70, 5.05, 5.30, 5.00, 4.75, 3.05, 7.05,
                6.00, 6.65, 5.85, 5.90, 6.90, 6.90, 3.85, 3.85,
                2.65, 4.45, 6.90, 5.30, 5.05, 4.05, 6.75, 5.00,
                5.50, 5.50, 9.10, 9.10, 6.65, 6.00, 5.90, 5.85),
        twist = (36.0, ) * 64,
        tilt =  (0.0, ) * 64,
        rise = 3.4,
    )
    # fmt: on

    name: str
    """Name of model."""

    order: int
    """Order of model, that is, length of oligonucleotides.

    Order 2 is a dinucleotide model, order 3 a trinucleotide model.
    """

    oligos: tuple[str]
    """All oligonucleotides of length :py:attr:`Model.order`."""

    rise: float
    """Displacement along the Z axis."""

    twist: dict[str, float]
    """Rotation angle in deg about the Z axis for all oligonucleotides."""

    roll: dict[str, float]
    """Rotation angle in deg about the Y axis for all oligonucleotides."""

    tilt: dict[str, float]
    """Rotation angle in deg about the X axis for all oligonucleotides."""

    matrices: dict[str | None, NDArray[Any]]
    """Homogeneous transformation matrices for all oligonucleotides."""

    def __init__(
        self,
        model: Model | dict[str, Any] | os.PathLike[Any] | str | None = None,
        /,
        **kwargs: Any,
    ) -> None:
        modeldict: dict[str, Any]
        if model:
            # TODO: type importfunction
            for importfunction in (
                self._fromname,
                self._fromdict,
                self._fromclass,
                self._fromfile,
            ):
                try:
                    # import functions return dictionary or raise exception
                    modeldict = importfunction(model)  # type: ignore
                    break
                except Exception:
                    pass
            else:
                raise ValueError(f'cannot initialize model from {model}')
        else:
            modeldict = Model.STRAIGHT

        modeldict.update(kwargs)

        try:
            self.oligos = modeldict['oligo'].split()
        except Exception:
            self.oligos = modeldict['oligo']
        self.order = len(self.oligos[0])
        self.name = str(modeldict['name'][:32])
        self.rise = float(modeldict['rise'])
        self.twist = dict(zip(self.oligos, modeldict['twist']))
        self.roll = dict(zip(self.oligos, modeldict['roll']))
        self.tilt = dict(zip(self.oligos, modeldict['tilt']))

        self.matrices = {}
        for oligo in oligonucleotides(self.order):
            if oligo not in self.twist:
                c = complementary(oligo)
                self.twist[oligo] = self.twist[c]
                self.roll[oligo] = self.roll[c]
                self.tilt[oligo] = -self.tilt[c]  # tilt reverses sign
            self.matrices[oligo] = dinucleotide_matrix(
                self.rise,
                self.twist[oligo],
                self.roll[oligo],
                self.tilt[oligo],
            ).T
        self.matrices[None] = dinucleotide_matrix(self.rise, 34.3, 0.0, 0.0).T

    def _fromfile(self, path: os.PathLike[Any] | str, /) -> dict[str, Any]:
        """Return model parameters as dict from file."""
        d: dict[str, Any] = {}
        with open(path, encoding='latin-1') as fh:
            d['name'] = fh.readline().rstrip()
            d['rise'] = float(fh.readline().split()[-1])

            def readtuple(
                itemtype: type, line: str, /
            ) -> tuple[tuple[Any, ...], str]:
                alist = [itemtype(i) for i in line.split()[1:]]
                while 1:
                    line = fh.readline()
                    if line.startswith('     '):
                        alist.extend(itemtype(i) for i in line.split())
                    else:
                        break
                return tuple(alist), line

            d['oligo'], line = readtuple(str, fh.readline())
            d['twist'], line = readtuple(float, line)
            d['roll'], line = readtuple(float, line)
            d['tilt'], line = readtuple(float, line)
        return d

    def _fromname(self, name: str, /) -> dict[str, Any]:
        """Return predefined model parameters as dict."""
        return getattr(Model, name.upper())  # type: ignore

    def _fromclass(self, aclass: Model, /) -> dict[str, Any]:
        """Return model parameters as dict from class."""
        return {a: getattr(aclass, a) for a in Model.STRAIGHT}

    def _fromdict(self, adict: dict[str, Any], /) -> dict[str, Any]:
        """Return model parameters as dict from dictionary."""
        for attr in Model.STRAIGHT:
            adict[attr]  # noqa: validation
        return adict

    def write(self, path: os.PathLike[Any] | str, /) -> None:
        """Write model to file."""
        with open(path, 'w', encoding='latin-1', newline='\n') as fh:
            fh.write(str(self))

    def save(self, path: os.PathLike[Any] | str, /) -> None:
        # deprecated
        self.write(path)

    def __repr__(self) -> str:
        return f'<{self.__class__.__name__} {self.name!r}>'

    def __str__(self) -> str:
        if self.order % 2:
            oligos = list(oligonucleotides(self.order))
        else:
            oligos = list(unique_oligos(self.order))

        def format_(
            items: Iterable[float | str],
            formatstr: str = '{:5.2f}',
            sep: str = '  ',
        ) -> str:
            i = [formatstr.format(item) for item in items]
            return '\n        '.join(sep.join(line) for line in chunks(i))

        return '\n'.join(
            (
                self.name.split('\n', maxsplit=1)[0],
                f'Rise    {self.rise:.2f}',
                'Oligo   ' + format_(oligos, '{}', ' ' * (7 - self.order)),
                'Twist   ' + format_(self.twist[i] for i in oligos),
                'Roll    ' + format_(self.roll[i] for i in oligos),
                'Tilt    ' + format_(self.tilt[i] for i in oligos),
            )
        )


class Sequence:
    r"""DNA nucleotide sequence.

    Parameters:
        arg: Sequence or name of file containing sequence.
        name: Name of sequence.
        comment: Single line description of sequence.
        maxlen: Maximum length of sequence.

    Notes:
        FASTA files must contain ``>name<space>comment<newline>sequence``.
        SEQ files must contain ``name<newline>comment<newline>sequence``.
        Nucleotides other than ATCG are ignored.

    Examples:
        >>> Sequence('0AxT-C:G a`t~c&g\t')[:]
        'ATCGATCG'
        >>> seq = Sequence('ATGCAAATTG' * 5, name='Test')
        >>> seq == 'ATGCAAATTG' * 5
        True
        >>> seq == None
        False
        >>> seq.write('_test.seq')
        >>> seq == Sequence('_test.seq')
        True

    """

    # Proc Natl Acad Sci USA, 1983, 80(24), p 7678, Fig 1
    KINETOPLAST = """
        GATCTAGACT AGACGCTATC GATAAAGTTT AAACAGTACA ACTATCGTGC TACTCACCTG
        TTGCCAAACA TTGCAAAAAT GCAAAATTGG GCTTGTGGAC GCGGAGAGAA TTCCCAAAAA
        TGTCAAAAAA TAGGCAAAAA ATGCCAAAAA TCCCAAACTT TTTAGGTCCC TCAGGTAGGG
        GCGTTCTCCG AAAACCGAAA AATGCATGCA GAAACCCCGT TCAAAAATCG GCCAAAATCG
        CCATTTTTTC AATTTTCGTG TGAAACTAGG GGTTGGTGTA AAATAGGGGT GGGGCTCCCC
        GGGGTAATTC TGGAAATTCG GGCCCTCAGG CTAGACCGGT CAAAATTAGG CCTCCTGACC
        CGTATATTTT TGGATTTCTA AATTTTGTGG CTTTAGATGT GGGAGATTTG GATC"""

    OUT_OF_PHASE_AAAAAA = 'CGCGCGCAAAAAACG'
    PHASED_AAAAAA = 'CGAAAAAACG'
    PHASED_GGGCCC = 'GAGGGCCCTA'

    name: str
    """Name of sequence."""

    comment: str
    """Single line description of sequence."""

    fname: str | None
    """File name."""

    _sequence: str

    def __init__(
        self,
        arg: os.PathLike[Any] | str,
        /,
        name: str = 'Untitled',
        comment: str = '',
        maxlen: int = 1024 * 1024,
    ) -> None:
        self.name = name if name else 'Untitled'
        self.comment = comment
        self._sequence = ''
        if os.path.isfile(arg):
            self._fromfile(arg, maxsize=maxlen * 2)
            self.fname = os.path.split(arg)[1]
        else:
            assert isinstance(arg, str)
            self._sequence = arg
            self.fname = None

        # clean name and comment
        for sep in (None, '|', ',', ';', '>', '<'):
            self.name = self.name.split(sep, 1)[0]
        self.name = self.name[:32]
        self.comment = (comment + ' ').splitlines()[0].strip()

        # remove all but ATCG from sequence
        nucls = dict(zip('ATCGatcg', 'ATCGATCG'))
        self._sequence = ''.join(nucls.get(c, '') for c in self._sequence)
        # limit length of sequence
        if maxlen and len(self._sequence) > maxlen:
            warnings.warn(f'sequence truncated to {maxlen} nucleotides')
            self._sequence = self._sequence[:maxlen]
        if not self._sequence:
            raise ValueError('not a valid sequence')

    def _fromfile(
        self, path: os.PathLike[Any] | str, /, maxsize: int = -1
    ) -> None:
        """Read name, comment and sequence from file."""
        with open(path, encoding='latin-1') as fh:
            firstline = fh.readline().rstrip()
            if firstline.startswith('>'):  # FASTA format
                self.name, self.comment = (firstline[1:] + ' ').split(' ', 1)
            elif firstline:
                self.name = firstline
                self.comment = fh.readline()
            self._sequence = fh.read(maxsize)

    def write(self, path: os.PathLike[Any] | str, /) -> None:
        """Write sequence to file.

        Parameters:
            path: Name of file to write.

        """
        with open(path, 'w', encoding='latin-1', newline='\n') as fh:
            fh.write(f'{self.name}\n{self.comment}\n{self.format()}')

    def save(self, path: os.PathLike[Any] | str, /) -> None:
        # deprecated
        self.write(path)

    @property
    def string(self) -> str:
        """Sequence string containing only ATCG."""
        return self._sequence

    def format(self, block: int = 10, line: int = 6) -> str:
        """Return string of sequence formatted in blocks and lines.

        Parameters:
            block: Length of blocks.
            Line: Line length.

        """
        blocks = chunks(self._sequence, block)
        lines = chunks(blocks, line)
        width = len(f'{(len(lines) - 1) * block * line}')
        return '\n'.join(
            '{:{}} {}'.format(i * line * block, width, ' '.join(s))
            for i, s in enumerate(lines)
        )

    def __getitem__(self, key: int | slice, /) -> str:
        return self._sequence[key]

    def __len__(self) -> int:
        return len(self._sequence)

    def __iter__(self) -> Iterator[str]:
        return iter(self._sequence)

    def __eq__(self, other: object, /) -> bool:
        if not isinstance(other, (str, Sequence)):
            return False
        return self._sequence == other[:]

    def __repr__(self) -> str:
        return (
            f'<{self.__class__.__name__} {self.name!r} '
            f'length={len(self._sequence)!r}>'
        )

    def __str__(self) -> str:
        return f'{self.name}\n{self.comment}\n{self.format()}'


def complementary(sequence: Sequence | str) -> str:
    """Return complementary DNA sequence.

    Parameters:
        sequence: DNA sequence.

    Examples:
        >>> complementary('AT CG')
        'CGAT'

    """
    c = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(sequence))


def oligonucleotides(
    length: int, /, nucleotides: str = 'AGCT'
) -> Iterator[str]:
    """Generate all oligonucleotide sequences of length.

    Parameters:
        length: Length of oligonucleotides to generate.
        nucleotides: Nucleotides in oligonucleotides.

    Yields:
        Oligonucletotide sequence of length.

    Examples:
        >>> ' '.join(oligonucleotides(2))
        'AA AG AC AT GA GG GC GT CA CG CC CT TA TG TC TT'

    """

    def rloop(length: int, part: str) -> Iterator[str]:
        if length:
            length -= 1
            for nucleotide in nucleotides:
                yield from rloop(length, part + nucleotide)
        else:
            yield part

    return rloop(length, '')


def unique_oligos(length: int, /, nucleotides: str = 'AGCT') -> Iterator[str]:
    """Generate all unique oligonucleotide sequences of length.

    Parameters:
        length: Length of oligonucleotides to generate.
        nucleotides: Nucleotides in oligonucleotides.

    Examples:
        >>> ' '.join(unique_oligos(2))
        'AA AG AC AT GA GG GC CA CG TA'

    """
    s: set[str] = set()
    for oligo in oligonucleotides(length, nucleotides):
        if oligo in s:
            s.remove(oligo)
        else:
            s.add(complementary(oligo))
            yield oligo


@overload
def chunks(sequence: str, size: int = 10, /) -> list[str]: ...


@overload
def chunks(sequence: list[str], size: int = 10, /) -> list[list[str]]: ...


def chunks(
    sequence: str | list[str],
    /,
    size: int = 10,
) -> list[str] | list[list[str]]:
    """Return sequence in chunks of size.

    Parameters:
        sequence: Sequence to be chunked.
        size: Length of chunks.

    Examples:
        >>> chunks('ATCG' * 4, 10)
        ['ATCGATCGAT', 'CGATCG']

    """
    return [  # type: ignore
        sequence[i : i + size] for i in range(0, len(sequence), size)
    ]


def overlapping_chunks(
    sequence: Sequence | str, size: int, overlap: int, /
) -> Iterator[tuple[int, str]]:
    """Return iterator over overlapping chunks of sequence.

    Parameters:
        sequence: Sequence to be chunked.
        size: Length of chunks.
        overlap: Size of overlap.

    Yields:
        Tuple of start position of chunk and chunk sequence.

    Examples:
        >>> list(overlapping_chunks('ATCG' * 4, 4, 2))
        [(0, 'ATCGATCG'), (4, 'ATCGATCG'), (8, 'ATCGATCG')]

    """
    index = 0
    while index < len(sequence) - 2 * overlap:
        yield index, sequence[index : index + size + 2 * overlap]
        index += size


def dinuc_window(
    sequence: Sequence | str, size: int, /
) -> Iterator[str | None]:
    """Return window of nucleotides around each dinucleotide in sequence.

    Parameters:
        sequence: Sequence to be windowed.
        size: Length of window.

    Yields:
        Oligonucleotide sequence at window or None if window overlaps border.

    Examples:
        >>> list(dinuc_window('ATCG', 2))
        ['AT', 'TC', 'CG']
        >>> list(dinuc_window('ATCG', 3))
        ['ATC', 'TCG', None]
        >>> list(dinuc_window('ATCG', 4))
        [None, 'ATCG', None]

    """
    assert 1 < size <= len(sequence)
    border = size // 2
    for i in range(0, border - 1):
        yield None
    for i in range(0, len(sequence) - size + 1):
        yield sequence[i : i + size]
    for i in range(len(sequence) - size + border, len(sequence) - 1):
        yield None


def dinucleotide_matrix(
    rise: float, twist: float, roll: float, tilt: float, /
) -> NDArray[Any]:
    """Return transformation matrix to move from one nucleotide to next.

    Parameters:
        rise: Displacement along the Z axis.
        twist: Rotation angle in deg about the Z axis.
        roll: Rotation angle in deg about the Y axis.
        tilt: Rotation angle in deg about the X axis.

    """
    twist = math.radians(twist)
    sinw = math.sin(twist)
    cosw = math.cos(twist)
    roll = math.radians(-roll)
    sinr = math.sin(roll)
    cosr = math.cos(roll)
    tilt = math.radians(-tilt)
    sint = math.sin(tilt)
    cost = math.cos(tilt)
    return numpy.array(
        (
            (
                cost * cosw,
                sinr * sint * cosw - cosr * sinw,
                cosr * sint * cosw + sinr * sinw,
                0.0,
            ),
            (
                cost * sinw,
                sinr * sint * sinw + cosr * cosw,
                cosr * sint * sinw - sinr * cosw,
                0.0,
            ),
            (-sint, sinr * cost, cosr * cost, rise),
            (0.0, 0.0, 0.0, 1.0),
        ),
        dtype=numpy.float64,
    )


def superimpose_matrix(v0: NDArray[Any], v1: NDArray[Any], /) -> NDArray[Any]:
    """Return matrix to transform given vector set to second vector set.

    Parameters:
        v0: Given vector set.
        v1: Target vector set.

    """
    # move centroids to origin
    t0 = numpy.mean(v0, axis=0)
    t1 = numpy.mean(v1, axis=0)
    v0 -= t0
    v1 -= t1
    # SVD of covariance matrix
    u, _, vh = numpy.linalg.svd(numpy.dot(v1.T, v0))
    # rotation matrix from SVD orthonormal bases
    R = numpy.dot(u, vh)
    if numpy.linalg.det(R) < 0.0:
        # correct reflections
        R -= numpy.outer(u[:, 2], vh[2, :] * 2.0)
    # homogeneous transformation matrix
    M = numpy.identity(4, dtype=numpy.float64)
    T = numpy.identity(4, dtype=numpy.float64)
    M[0:3, 0:3] = R
    M[:3, 3] = t1
    T[0:3, 3] = -t0
    return numpy.dot(M, T)


def norm(vector: NDArray[Any], /) -> float:
    """Return length of vector, that is, its euclidean norm.

    Parameters:
        vector: Vector.

    """
    # return numpy.linalg.norm(vector)
    return float(numpy.sqrt(numpy.dot(vector, vector)))


def main(argv: list[str] | None = None, /) -> int:
    """Command line usage main function.

    Parameters:
        argv: Command line arguments.

    """
    if argv is None:
        argv = sys.argv

    # TODO: use argparse module
    import optparse

    def search_doc(r: str, d: str) -> str:
        if not __doc__:
            return d
        match = re.search(r, __doc__)
        if match is None:
            return d
        return match.group(1)

    parser = optparse.OptionParser(
        usage='usage: %prog [options] sequence | file',
        description=search_doc('\n\n([^|]*?)\n\n', ''),
        version=f'%prog {__version__}',
        prog='dnacurve',
    )
    opt = parser.add_option
    opt(
        '-m',
        '--model',
        dest='model',
        metavar='MODEL',
        default='trifonov',
        help='input model name or file',
    )
    opt(
        '-n',
        '--name',
        dest='name',
        metavar='NAME',
        default='Untitled',
        help='set sequence name',
    )
    opt('--csv', dest='csv', metavar='FILE', help='save results as CSV')
    opt('--pdb', dest='pdb', metavar='FILE', help='save coordinates as PDB')
    opt('--seq', dest='seq', metavar='FILE', help='save nucleotide sequence')
    opt('--png', dest='png', metavar='FILE', help='save plot as PNG image')
    opt('--pdf', dest='pdf', metavar='FILE', help='save plot as PDF')
    opt('--ps', dest='ps', metavar='FILE', help='save plot as PostScript')
    opt(
        '-p',
        '--plot',
        dest='plot',
        action='store_true',
        default=False,
        help='plot to interactive window',
    )
    opt(
        "--dpi", dest='dpi', type='int', default=96, help='set plot resolution'
    )
    opt(
        '--web',
        dest='web',
        action='store_true',
        default=False,
        help='start web application and open it in a web browser',
    )
    opt(
        '--url',
        dest='url',
        default=None,
        help='URL to run web application',
    )
    opt(
        '--nobrowser',
        dest='nobrowser',
        action='store_true',
        default=False,
        help='do not open web browser',
    )
    opt(
        '--test',
        dest='test',
        action='store_true',
        default=False,
        help='analyze a test sequence',
    )
    opt(
        '--doctest',
        dest='doctest',
        action='store_true',
        default=False,
        help='run the internal tests',
    )
    opt('-v', '--verbose', dest='verbose', action='store_true', default=True)
    opt('-q', '--quiet', dest='verbose', action='store_false')

    settings, sequence_list = parser.parse_args()

    if settings.web:
        try:
            from .web import main as web_main
        except ImportError:
            from web import main as web_main  # type: ignore

        return web_main(settings.url, not settings.nobrowser)
    if settings.doctest:
        import doctest

        numpy.set_printoptions(suppress=True, precision=5)
        try:
            import dnacurve.dnacurve as m
        except ImportError:
            m = None  # type: ignore
        doctest.testmod(m, optionflags=doctest.ELLIPSIS)
        return 0
    if settings.test:
        settings.name = 'Kinetoplast'
        sequence = Sequence.KINETOPLAST
    elif sequence_list:
        settings.name = settings.name[:16]
        sequence = ''.join(sequence_list)
    else:
        parser.error('no sequence specified')

    def file_ext(f: str | None, e: str) -> str | None:
        return f if (f is None) else (f if f.lower().endswith(e) else f + e)

    settings.csv = file_ext(settings.csv, '.csv')
    settings.seq = file_ext(settings.seq, '.seq')
    settings.pdb = file_ext(settings.pdb, '.pdb')
    settings.png = file_ext(settings.png, '.png')
    settings.pdf = file_ext(settings.pdf, '.pdf')
    settings.ps = file_ext(settings.ps, '.ps')

    try:
        results = CurvedDNA(sequence, settings.model, name=settings.name)
    except Exception as exc:
        print('\nFatal Error: \n  ', exc, sep='')
        # raise
    else:
        if settings.verbose:
            print('\n', results, sep='')
        if settings.csv:
            results.write_csv(settings.csv)
        if settings.pdb:
            results.write_pdb(settings.pdb)
        if settings.seq:
            results.sequence.write(settings.seq)
        if settings.plot:
            try:
                results.plot(settings.plot, dpi=settings.dpi)
            except Exception as exc:
                print('\nFailed to plot results: \n  ', exc, sep='')
        else:
            if settings.ps:
                results.plot(settings.ps, dpi=settings.dpi)
            if settings.pdf:
                results.plot(settings.pdf, dpi=settings.dpi)
            if settings.png:
                results.plot(settings.png, dpi=settings.dpi)
    return 0


MODELS.extend(
    sorted(
        (a for a in dir(Model) if not a.startswith('_') and a.isupper()),
        key=lambda x: getattr(Model, x)['name'],  # type: ignore
    )
)

if __name__ == '__main__':
    sys.exit(main())
