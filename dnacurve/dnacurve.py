# dnacurve.py

# Copyright (c) 1994-2020, Christoph Gohlke
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

"""DNA Curvature Analysis.

Dnacurve is a Python library and console script to calculate the global
3D structure of a B-DNA molecule from its nucleotide sequence according to the
dinucleotide wedge model. Local bending angles and macroscopic curvature
are calculated at each nucleotide.

For command line usage run ``python -m dnacurve --help``

:Author: `Christoph Gohlke <https://www.lfd.uci.edu/~gohlke/>`_

:License: BSD 3-Clause

:Version: 2020.1.1

Requirements
------------
* `CPython >= 3.6 <https://www.python.org>`_
* `Numpy 1.14 <https://www.numpy.org>`_
* `Matplotlib 3.1 <https://www.matplotlib.org>`_

Revisions
---------
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
The API is not stable yet and is expected to change between revisions.

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
>>> result = CurvedDNA('ATGCAAATTG'*5, 'trifonov', name='Example')
>>> result.curvature[:, 18:22]
array([[0.58062, 0.58163, 0.58278, 0.58378],
       [0.0803 , 0.11293, 0.07676, 0.03166],
       [0.57924, 0.5758 , 0.57368, 0.5735 ]])
>>> result.save_csv('_test.csv')
>>> result.save_pdb('_test.pdb')
>>> result.plot('_test.png', dpi=160)

"""

__version__ = '2020.1.1'

__all__ = (
    'CurvedDNA', 'Model', 'Sequence', 'MODELS', 'MAXLEN', 'main',
    'complementary', 'oligonucleotides', 'unique_oligos', 'chunks',
    'overlapping_chunks', 'dinuc_window', 'dinucleotide_matrix',
    'superimpose_matrix'
)

import sys
import os
import re
import math
import datetime
import warnings

import numpy

MAXLEN = 510  # maximum length of sequence


class CurvedDNA:
    """Calculate, plot or save helix coordinates, local bending and curvature.

    Attributes
    ----------
    model : Model
        The curvature model.
    sequence : Sequence
        The DNA sequence to analyze.
    coordinates : 3D ndarray
        Homogeneous coordinates at each nucleotide of:
        Index 0) helix axis,
        Index 1) phosphate of 5'-3' strand,
        Index 2) phosphate of antiparallel strand,
        Index 3) basepair normal vector,
        Index 4) smoothed basepair normal vector.
    curvature : 2D ndarray
        Values at each nucleotide, normalized relative to curvature in
        nucleosome:
        Index 0) curvature,
        Index 1) local bend angle,
        Index 2) curvature angle.
    windows : sequence of int
        Window sizes for calculating curvature, local bend angle, and
        curvature angle.
    scales : 2D ndarray
        Scaling factors used to normalize curvature array.

    Notes
    -----
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

    Examples
    --------
    See module examples.

    """

    P_COORDINATES = (  # cylindrical coordinates of 5' phosphate
        8.91,  # distance from axis
        -5.2,  # angle to roll axis
        2.08,  # distance from basepair plane
    )

    def __init__(self, sequence, model='trifonov', name='Untitled',
                 curvature_window=10, bend_window=2, curve_window=15,
                 maxlen=MAXLEN):
        """Initialize instance from sequence and model.

        Parameters
        ----------
        sequence : various types
            Sequence instance, file name, or nucleotide sequence.
            See Sequence constructor documentation.
        model : various types
            Model instance, file name, class, dict, or name of
            predefined model. See Model constructor documentation.
        name: str
            Optional human readable label.
        curvature_window : int
            Window size for calculating the curvature (default 10).
        bend_window : int
            Window size for calculating local bend angles (default 2).
        curve_window : int
            Window size for calculating curvature angles (default 15).
        maxlen : int
            Maximum length of sequence (default: 500)

        """
        if isinstance(model, Model):
            self.model = model
        else:
            self.model = Model(model)
        if isinstance(sequence, Sequence):
            self.sequence = sequence
        else:
            self.sequence = Sequence(sequence, name, maxlen=maxlen)
        if len(self.sequence) < self.model.order:
            raise ValueError(
                f'sequence must be >{self.model.order} nucleotides long'
            )
        if len(self.sequence) > maxlen:
            warnings.warn(f'sequence is longer than {maxlen} nucleotides')

        assert 0 < curvature_window < 21
        assert 0 < bend_window < 4
        assert 9 < curve_window < 21
        self.windows = [curvature_window, bend_window, curve_window]
        self._limits = [10.0, 10.0, 10.0]

        self.date = datetime.datetime.now()
        self.coordinates = numpy.zeros((5, len(self), 4), dtype='float64')
        self.curvature = numpy.zeros((3, len(self)), dtype='float64')
        self.scales = numpy.ones((3, 1), dtype='float64')

        self._coordinates()
        self._reorient()
        self._center()
        self._curvature()

    def __len__(self):
        """Return number of nucleotides in sequence."""
        return len(self.sequence)

    def __str__(self):
        """Return string representation of sequence and model."""
        return f'{self.sequence}\n\n{self.model}\n'

    def _coordinates(self):
        """Calculate coordinates and normal vectors from sequence and model."""
        p = self.P_COORDINATES
        p = numpy.array((p[0] * math.cos(math.radians(p[1])),
                         p[0] * math.sin(math.radians(p[1])),
                         p[2]))

        xyz = self.coordinates
        xyz[0:3, :, 3] = 1.0  # homogeneous coordinates
        xyz[1, :, 0:3] = p  # 5' phosphate
        xyz[2, :, 0:3] = -p[0], p[1], -p[2]  # phosphate of antiparallel strand
        xyz[3, :, 2] = 1.0  # basepair normal vectors

        matrices = self.model.matrices
        dot = numpy.dot
        for i, seq in enumerate(dinuc_window(self.sequence, self.model.order)):
            xyz[:4, :i + 1, :] = dot(xyz[:4, :i + 1, :], matrices[seq])

        # Average direction vector of one helix turn,
        # calculated by smoothing the basepair normals
        if len(self.sequence) > 10:
            kernel = numpy.array([0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5])
            kernel /= kernel.sum()
            for i in 0, 1, 2:
                xyz[4, :, i] = numpy.convolve(kernel, xyz[3, :, i], 'same')
            for i in range(5, len(self) - 5):
                xyz[4, i, :] /= norm(xyz[4, i, :])

    def _reorient(self):
        """Reorient coordinates."""
        xyz = self.coordinates[0, :, 0:3]  # helix axis
        xyz = xyz - xyz[-1]
        # assert start point is at origin
        assert numpy.allclose(xyz[-1], (0, 0, 0))
        # normalized end to end vector
        e = +xyz[0]
        e_len = norm(e)
        e /= e_len
        # point i of maximum distance to end to end line
        x = numpy.cross(e, xyz)
        x = numpy.sum(x * x, axis=1)
        i = numpy.argmax(x)
        x = math.sqrt(x[i])
        # distance of endpoint to xyz[i]
        w = norm(xyz[i])
        # distance of endpoint to point on end to end line nearest to xyz[i]
        u = math.sqrt(w * w - x * x)
        # find transformation matrix
        v0 = xyz[[0, i, -1]]
        v1 = numpy.array(((0, 0, 0), (e_len - u, 0, x), (e_len, 0, 0)))
        M = superimpose_matrix(v0, v1)
        self.coordinates = numpy.dot(self.coordinates, M.T)

    def _center(self):
        """Center atomic coordinates at origin."""
        xyz = self.coordinates[0:3, :, 0:3]  # helix axis and P atoms
        low = numpy.min(numpy.min(xyz, axis=1), axis=0)
        upp = numpy.max(numpy.max(xyz, axis=1), axis=0)
        self._limits = (upp - low) / 2.0
        self.coordinates[0:3, :, 0:3] -= low + self._limits

    def _curvature(self):
        """Calculate normalized curvature and bend angles."""
        dot = numpy.dot
        cross = numpy.cross
        arccos = numpy.arccos

        # curvature from radius
        window = self.windows[0]
        if len(self) >= 2 * window:
            result = self.curvature[0, :]
            xyz = self.coordinates[0, :, 0:3]  # helix axis
            for i in range(window, len(self) - window):
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
        if len(self) >= 2 * window:
            normals = self.coordinates[3, :, 0:3]
            result = self.curvature[1, :]
            for i in range(window, len(self) - window):
                if not numpy.allclose(normals[i - window],
                                      normals[i + window]):
                    result[i] = arccos(dot(normals[i - window],
                                           normals[i + window]))

        # curvature angles from smoothed basepair normals
        window = self.windows[2]
        if len(self) >= 2 * window:
            normals = self.coordinates[4, :, 0:3]
            result = self.curvature[2, :]
            for i in range(window, len(self) - window):
                if not numpy.allclose(normals[i - window],
                                      normals[i + window]):
                    result[i] = arccos(dot(normals[i - window],
                                           normals[i + window]))

        self.curvature = numpy.nan_to_num(self.curvature)

        # normalize relative to curvature in nucleosome
        self.scales[0] = 0.0234
        self.scales[1:] = 0.0234 * 2 * self.windows[2] * self.model.rise
        self.curvature /= self.scales

    def save_csv(self, path):
        """Save coordinates and curvature values to CSV file."""
        with open(path, 'w') as fh:
            fh.write(self.csv())

    def save_pdb(self, path):
        """Save atomic coordinates to PDB file."""
        with open(path, 'w') as fh:
            fh.write(self.pdb())

    def csv(self):
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
                self.scales[0, 0], self.windows[0],
                self.scales[1, 0], self.windows[1],
                self.scales[2, 0], self.windows[2],
            )
        ]
        for i in range(len(self)):
            csv.append(f'{seq[i]},{i + 1}')
            csv.append(',{:.5f},{:.5f},{:.5f}'.format(*cur[:, i]))
            for j in range(5):
                csv.append(',{:.5f},{:.5f},{:.5f}'.format(*xyz[j, i, 0:3]))
            csv.append('\n')
        return ''.join(csv)

    def pdb(self):
        """Return atomic coordinates in PDB format."""
        pdb = []
        date1 = datetime.datetime.strftime(self.date, '%d-%b-%y').upper()
        date2 = datetime.datetime.strftime(self.date, '%d %b %Y, %H:%M:%S')
        pdb.append(f'HEADER    {self.sequence.name:<40s}{date1:9s}   0XXX\n')
        pdb.append(f'REMARK    Generated by dnacurve.py on {date2}\n')
        pdb.append(f'REMARK    Sequence: {self.sequence.name}\n')
        pdb.append(f'REMARK    Parameters: {self.model.name}\n')
        seq = self.sequence
        xyz = self.coordinates[0:3, :, 0:3]
        xyz = xyz - xyz.min(axis=1).min(axis=0)
        for j, (name, chain, element) in enumerate((('XA', 'A', 'C'),
                                                    ('P', 'B', 'P'),
                                                    ('P', 'C', 'P'))):
            for i in range(len(self)):
                x, y, z = xyz[j, i, 0:3]
                pdb.append(
                    f'ATOM  {i + 1:5} {name:<4s} {seq[i]:<3s}'
                    f' {chain:1s}{i + 1:4}    {x:8.3f}{y:8.3f}{z:8.3f}'
                    f'  1.00  0.00          {element:>2s}  \n'
                )
            pdb.append('TER\n')
        pdb.append('END\n')
        return ''.join(pdb)

    def plot(self, arg=True, dpi=96, figsize=(6.0, 7.5), imageformat=None):
        """Plot results using matplotlib.

        Parameters
        ----------
        arg : bool or str
            False: do not plot.
            True: interactive plot.
            String: path name to save figure.
        dpi : int
            Resolution of plot in dots per inch.

        """
        if not arg:
            return None

        if arg is True:
            # GUI plot
            from matplotlib import pyplot

            fig = pyplot.figure(dpi=dpi, figsize=figsize)
            try:
                fig.canvas.manager.window.title('DNA Curvature Analysis')
            except AttributeError:
                pass
        else:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            from matplotlib.figure import Figure

            fig = Figure(dpi=dpi, figsize=figsize)
            FigureCanvasAgg(fig)

        # sequence
        ax = fig.add_subplot(321)
        ax.set_title(
            f'Sequence: {self.sequence.name}\nModel: {self.model.name}',
            x=0, horizontalalignment='left', size=11
        )
        ax.text(
            0., 0.95, self.sequence.format(line=3), family='monospace',
            size=7, verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes
        )
        ax.set_axis_off()
        ax.axis('image')

        # projections
        xyz = self.coordinates[..., 0:3]
        limit = numpy.max(self._limits) * 1.1

        def plot_projection(plotnum, axes, label=True):
            ax = fig.add_subplot(plotnum)
            ax.set_title(f'{axes[0]}-{axes[1]}', size=11)
            ax0, ax1 = (ord(a) - 88 for a in axes)
            # TODO: sort lines by depth and draw back to front
            ax.plot(xyz[1, :, ax0], xyz[1, :, ax1], 'r-',
                    xyz[2, :, ax0], xyz[2, :, ax1], 'b-',
                    xyz[0, :, ax0], xyz[0, :, ax1], 'k-', lw=0.8)
            if label:
                ax.text(xyz[1, 0, ax0], xyz[1, 0, ax1], "5'",
                        color='darkred', clip_on=True, size=10)
                ax.text(xyz[1, -1, ax0], xyz[1, -1, ax1], "3'",
                        color='darkred', clip_on=True, size=10)
            ax.axis('image')
            ax.axis([-limit, limit, -limit, limit])
            ax.axis('off')

        plot_projection(322, 'XZ')
        plot_projection(323, 'ZY', False)
        plot_projection(324, 'XY')

        # plot
        ax = fig.add_subplot(325)
        ax.set_position((0.1, 0.075, 0.83, 0.21))

        def plot_curvature(index, label, style, lw):
            w = self.windows[index]
            ax.plot(numpy.arange(w, len(self) - w),
                    self.curvature[index, w:-w], style, lw=lw,
                    label=label + f' / {self.scales[index, 0]:.4f} [{w}]')

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
        ax.axis((0, len(self), 0, 1.0))
        ax.set_xlabel('Sequence Index', size=11)

        if arg is True:
            pyplot.show()
        elif arg:
            fig.savefig(arg, dpi=dpi, format=imageformat)
        else:
            return fig
        return None

    @property
    def name(self):
        """Return name of sequence."""
        return self.sequence.name


class Model:
    """N-mer DNA-bending model.

    Transformation parameters and matrices for all oligonucleotides of
    certain length.

    Attributes
    ----------
    name : str
        Human readable label.
    order : int
        Order of model, i.e. length of oligonucleotides.
        Order 2 is a dinucleotide model, order 3 a trinucleotide model etc.
    rise : float
        Displacement along the Z axis.
    twist : dict
        Rotation angle in deg about the Z axis for all oligonucleotides.
    roll : dict
        Rotation angle in deg about the Y axis for all oligonucleotides.
    tilt : dict
        Rotation angle in deg about the Z axis for all oligonucleotides.
    matrices : dict
        Homogeneous transformation matrices for all oligonucleotides.

    Examples
    --------
    >>> m = Model('AAWedge')
    >>> m = Model('Nucleosome')
    >>> m = Model(**Model.STRAIGHT)
    >>> m = Model(Model.CALLADINE, name='My Model', rise=4.0)
    >>> m.name == 'My Model' and m.rise == 4.0
    True
    >>> m = Model(name='Test', rise=3.38,
    ...           oligo='AA AC AG AT CA GG CG GA GC TA'.split(),
    ...           twist=(34.29, )*10, roll=(0., )*10, tilt=(0., )*10)
    >>> m.save('_test.dat')
    >>> m.twist == Model('_test.dat').twist
    True

    """

    STRAIGHT = dict(
        name = 'Straight',
        oligo = 'AA AC AG AT CA GG CG GA GC TA',
        twist = (360.0 / 10.5, ) * 10,
        roll =  (0.0, ) * 10,
        tilt =  (0.0, ) * 10,
        rise = 3.38)

    AAWEDGE = dict(
        name = 'AA Wedge',
        oligo = 'AA     AC    AG    AT    CA    GG     CG    GA    GC    TA',
        twist = (35.62, 34.4, 27.7, 31.5, 34.5, 33.67, 29.8, 36.9, 40.0, 36.0),
        roll =  (-8.40,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0),
        tilt =  ( 2.40,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0),
        rise = 3.38)

    CALLADINE = dict(
        # Nucleic Acids Res, 1994, 22(24), p 5498, Table 1, Model b.
        name = 'Calladine & Drew',
        oligo = 'AA    AC    AG    AT    CA    GG    CG    GA    GC    TA',
        twist = (35.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0),
        roll =  ( 0.0,  3.3,  3.3,  3.3,  3.3,  3.3,  3.3,  3.3,  3.3,  6.6),
        tilt =  ( 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0),
        rise = 3.38)

    TRIFONOV = dict(
        # Nucleic Acids Res, 1994, 22(24), p 5498, Table 1, Model c.
        name = 'Bolshoi & Trifonov',
        oligo = 'AA     AC    AG    AT    CA    GG     CG    GA    GC    TA',
        twist = (35.62, 34.4, 27.7, 31.5, 34.5, 33.67, 29.8, 36.9, 40.0, 36.0),
        roll =  (-6.50, -0.9,  8.4,  2.6,  1.6,   1.2,  6.7, -2.7, -5.0,  0.9),
        tilt =  ( 3.20, -0.7, -0.3,  0.0,  3.1,  -1.8,  0.0, -4.6,  0.0,  0.0),
        rise = 3.38)

    DESANTIS = dict(
        # Nucleic Acids Res, 1994, 22(24), p 5498, Table 1, Model d.
        name = 'Cacchione & De Santis',
        oligo = 'AA    AC    AG    AT    CA    GG    CG    GA    GC    TA',
        twist = (35.9, 34.6, 35.6, 35.0, 34.5, 33.0, 33.7, 35.8, 33.3, 34.6),
        roll =  (-5.4, -2.4,  1.0, -7.3,  6.7,  1.3,  4.6,  2.0, -3.7,  8.0),
        tilt =  (-0.5, -2.7, -1.6,  0.0,  0.4, -0.6,  0.0, -1.7,  0.0,  0.0),
        rise = 3.38)

    REVERSED = dict(
        # Nucleic Acids Res, 1994, 22(24), p 5498, Table 1, Model e.
        name = 'Reversed Calladine & Drew',
        oligo = 'AA    AC    AG    AT    CA    GG    CG    GA    GC    TA',
        twist = (35.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0, 34.0),
        roll =  ( 3.3,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -3.3),
        tilt =  ( 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0),
        rise = 3.38)

    NUCLEOSOME = dict(
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
        rise = 3.38)

    TRINUCLEOTIDE = dict(
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
        rise = 3.4)

    def __init__(self, model=None, **kwargs):
        """Initialize instance from predefined model, file, or arguments.

        Parameters
        ----------
        model : various types
            Name of predefined model : str
                'straight', 'aawedge', 'trifonov', 'desantis', 'calladine',
                'reversed'
            Class or Dict:
                Instance containing model parameters
            Path name: str
                File containing model parameters
            None:
                Default model 'straight'
        name : str
            Human readable label.
        oligo : str or tuple
            Oligonucleotide sequences separated by whitespace or as tuple.
        twist : sequence of floats
            Twist values for given oligonucleotides in degrees.
        roll : sequence of floats
            Roll values for given oligonucleotides in degrees.
        tilt : sequence of floats
            Tilt values for given oligonucleotides in degrees.
        rise : float
            Rise value.

        """
        if model:
            for importfunction in (self._fromname,
                                   self._fromdict,
                                   self._fromclass,
                                   self._fromfile):
                try:
                    # import functions return dictionary or raise exception
                    model = importfunction(model)
                    break
                except Exception:
                    pass
            else:
                raise ValueError(f'cannot initialize model from {model}')
        else:
            model = Model.STRAIGHT

        model.update(kwargs)

        try:
            self.oligos = model['oligo'].split()
        except Exception:
            self.oligos = model['oligo']
        self.order = len(self.oligos[0])
        self.name = str(model['name'][:32])
        self.rise = float(model['rise'])
        self.twist = dict(zip(self.oligos, model['twist']))
        self.roll = dict(zip(self.oligos, model['roll']))
        self.tilt = dict(zip(self.oligos, model['tilt']))

        self.matrices = {}
        for oligo in oligonucleotides(self.order):
            if oligo not in self.twist:
                c = complementary(oligo)
                self.twist[oligo] = self.twist[c]
                self.roll[oligo] = self.roll[c]
                self.tilt[oligo] = -self.tilt[c]  # tilt reverses sign
            self.matrices[oligo] = dinucleotide_matrix(
                self.rise, self.twist[oligo], self.roll[oligo],
                self.tilt[oligo]).T
        self.matrices[None] = dinucleotide_matrix(self.rise, 34.3, 0.0, 0.0).T

    def __str__(self):
        """Return string representation of model."""
        if self.order % 2:
            oligos = list(oligonucleotides(self.order))
        else:
            oligos = list(unique_oligos(self.order))

        def format_(items, formatstr='{:5.2f}', sep='  '):
            items = [formatstr.format(item) for item in items]
            return '\n        '.join(sep.join(line) for line in chunks(items))

        return '\n'.join((
            '{}'.format(self.name.split('\n')[0]),
            'Rise    {:.2f}'.format(self.rise),
            'Oligo   ' + format_(oligos, '{}', ' ' * (7 - self.order)),
            'Twist   ' + format_(self.twist[i] for i in oligos),
            'Roll    ' + format_(self.roll[i] for i in oligos),
            'Tilt    ' + format_(self.tilt[i] for i in oligos),
        ))

    def _fromfile(self, path):
        """Return model parameters as dict from file."""
        d = {}
        with open(path, 'r') as fh:
            d['name'] = fh.readline().rstrip()
            d['rise'] = float(fh.readline().split()[-1])

            def readtuple(itemtype, line):
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

    def _fromname(self, name):
        """Return predefined model parameters as dict."""
        return getattr(Model, name.upper())

    def _fromclass(self, aclass):
        """Return model parameters as dict from class."""
        return {a: getattr(aclass, a) for a in Model.STRAIGHT}

    def _fromdict(self, adict):
        """Return model parameters as dict from dictionary."""
        for attr in Model.STRAIGHT:
            adict[attr]  # noqa: validation
        return adict

    def save(self, path):
        """Save model to file."""
        with open(path, 'w') as fh:
            fh.write(str(self))


class Sequence:
    """DNA nucleotide sequence.

    Attributes
    ----------
    name : str
        Human readable label.
    comment : str
        Single line description of sequence.
    string : str
        Sequence string containing only ATCG.

    Notes
    -----
    FASTA files must contain ``>name<space>comment<newline>sequence``.
    SEQ files must contain ``name<newline>comment<newline>sequence``.
    Nucleotides other than ATCG are ignored.

    Examples
    --------
    >>> Sequence('0AxT-C:G a`t~c&g\t')[:]
    'ATCGATCG'
    >>> seq = Sequence('ATGCAAATTG'*3, name='Test')
    >>> seq == 'ATGCAAATTG'*3
    True
    >>> seq == None
    False
    >>> seq.save('_test.seq')
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

    def __init__(self, arg, name='Untitled', comment='', maxlen=1024*1024):
        """Initialize instance from nucleotide sequence string or file name."""
        self.name = name
        self.comment = comment
        self._sequence = ''
        if os.path.isfile(arg):
            self._fromfile(arg, maxsize=maxlen * 2)
            self.fname = os.path.split(arg)[1]
        else:
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

    def _fromfile(self, path, maxsize=-1):
        """Read name, comment and sequence from file."""
        with open(path, 'r') as fh:
            firstline = fh.readline().rstrip()
            if firstline.startswith('>'):  # Fasta format
                self.name, self.comment = (firstline[1:] + ' ').split(' ', 1)
            elif firstline:
                self.name = firstline
                self.comment = fh.readline()
            self._sequence = fh.read(maxsize)

    def save(self, path):
        """Save sequence to file."""
        with open(path, 'w') as fh:
            fh.write(str(self))

    @property
    def string(self):
        """Return sequence as string."""
        return self._sequence

    def format(self, block=10, line=6):
        """Return string of sequence formatted in blocks and lines."""
        lines = chunks(chunks(self._sequence, block), line)
        width = len(f'{(len(lines) - 1) * block * line}')
        for i, s in enumerate(lines):
            lines[i] = '{:{}} {}'.format(i * line * block, width, ' '.join(s))
        return '\n'.join(lines)

    def __getitem__(self, key):
        """Return nucleotide at position."""
        return self._sequence[key]

    def __len__(self):
        """Return number of nucleotides in the sequence."""
        return len(self._sequence)

    def __iter__(self):
        """Return iterator over nucleotides."""
        return iter(self._sequence)

    def __eq__(self, other):
        """Return result of sequence comparison."""
        try:
            return self._sequence == other[:]
        except Exception:
            return False

    def __str__(self):
        """Return string representation of sequence."""
        return f'{self.name}\n{self.comment}\n{self.format()}'


def complementary(sequence):
    """Return complementary sequence.

    Examples
    --------
    >>> complementary('AT CG')
    'CGAT'

    """
    c = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(sequence))


def oligonucleotides(length, nucleotides='AGCT'):
    """Generate all oligonucleotide sequences of length.

    Examples
    --------
    >>> ' '.join(oligonucleotides(2))
    'AA AG AC AT GA GG GC GT CA CG CC CT TA TG TC TT'

    """

    def rloop(length, part):
        if length:
            length -= 1
            for nucleotide in nucleotides:
                yield from rloop(length, part + nucleotide)
        else:
            yield part

    return rloop(length, '')


def unique_oligos(length, nucleotides='AGCT'):
    """Generate all unique oligonucleotide sequences of length.

    Examples
    --------
    >>> ' '.join(unique_oligos(2))
    'AA AG AC AT GA GG GC CA CG TA'

    """
    s = set()
    for oligo in oligonucleotides(length, nucleotides):
        if oligo in s:
            s.remove(oligo)
        else:
            s.add(complementary(oligo))
            yield oligo


def chunks(sequence, size=10):
    """Return sequence in chunks of size.

    Examples
    --------
    >>> chunks('ATCG'*4, 10)
    ['ATCGATCGAT', 'CGATCG']

    """
    return [sequence[i:i + size] for i in range(0, len(sequence), size)]


def overlapping_chunks(sequence, size, overlap):
    """Return iterator over overlapping chunks of sequence.

    Examples
    --------
    >>> list(overlapping_chunks('ATCG'*4, 4, 2))
    [(0, 'ATCGATCG'), (4, 'ATCGATCG'), (8, 'ATCGATCG')]

    """
    index = 0
    while index < len(sequence) - 2 * overlap:
        yield index, sequence[index:index+size+2*overlap]
        index += size


def dinuc_window(sequence, size):
    """Return window of nucleotides around each dinucleotide in sequence.

    Return None if window overlaps border.

    Examples
    --------
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
        yield sequence[i:i + size]
    for i in range(len(sequence) - size + border, len(sequence) - 1):
        yield None


def dinucleotide_matrix(rise, twist, roll, tilt):
    """Return transformation matrix to move from one nucleotide to next."""
    twist = math.radians(twist)
    sinw = math.sin(twist)
    cosw = math.cos(twist)
    roll = math.radians(-roll)
    sinr = math.sin(roll)
    cosr = math.cos(roll)
    tilt = math.radians(-tilt)
    sint = math.sin(tilt)
    cost = math.cos(tilt)
    return numpy.array((
        (cost*cosw, sinr*sint*cosw-cosr*sinw, cosr*sint*cosw+sinr*sinw,  0.0),
        (cost*sinw, sinr*sint*sinw+cosr*cosw, cosr*sint*sinw-sinr*cosw,  0.0),
        (    -sint,                sinr*cost,                cosr*cost, rise),
        (      0.0,                      0.0,                      0.0,  1.0)),
                       dtype='float64')


def superimpose_matrix(v0, v1):
    """Return matrix to transform given vector set to second vector set."""
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
    M = numpy.identity(4, dtype='float64')
    T = numpy.identity(4, dtype='float64')
    M[0:3, 0:3] = R
    M[:3, 3] = t1
    T[0:3, 3] = -t0
    return numpy.dot(M, T)


def norm(vector):
    """Return length of vector, i.e. its euclidean norm."""
    # return numpy.linalg.norm(vector)
    return numpy.sqrt(numpy.dot(vector, vector))


def main(argv=None):
    """Command line usage main function."""
    if argv is None:
        argv = sys.argv

    # TODO: use argparse module
    import optparse

    def search_doc(r, d):
        return re.search(r, __doc__).group(1) if __doc__ else d

    parser = optparse.OptionParser(
        usage='usage: %prog [options] sequence | file',
        description=search_doc('\n\n([^|]*?)\n\n', ''),
        version=f'%prog {__version__}',
        prog='dnacurve',
    )
    opt = parser.add_option
    opt('-m', '--model', dest='model', metavar='MODEL', default='trifonov',
        help='input model name or file')
    opt('-n', '--name', dest='name', metavar='NAME', default='Untitled',
        help='set sequence name')
    opt('--csv', dest='csv', metavar='FILE', help='save results as CSV')
    opt('--pdb', dest='pdb', metavar='FILE', help='save coordinates as PDB')
    opt('--seq', dest='seq', metavar='FILE', help='save nucleotide sequence')
    opt('--png', dest='png', metavar='FILE', help='save plot as PNG image')
    opt('--pdf', dest='pdf', metavar='FILE', help='save plot as PDF')
    opt('--ps', dest='ps', metavar='FILE', help='save plot as PostScript')
    opt('-p', '--plot', dest='plot', action='store_true', default=False,
        help='plot to interactive window')
    opt("--dpi", dest='dpi', type='int', default=96,
        help='set plot resolution')
    opt('--web', dest='web', action='store_true', default=False,
        help='start web application and open it in a web browser')
    opt('--test', dest='test', action='store_true', default=False,
        help='analyze a test sequence')
    opt('--doctest', dest='doctest', action='store_true', default=False,
        help='run the internal tests')
    opt('-v', '--verbose', dest='verbose', action='store_true', default=True)
    opt('-q', '--quiet', dest='verbose', action='store_false')

    settings, sequence = parser.parse_args()

    if settings.web:
        try:
            from . import dnacurve_web  # noqa: delay import
        except ImportError:
            import dnacurve_web  # noqa: delay import

        return dnacurve_web.main()
    if settings.doctest:
        import doctest

        numpy.set_printoptions(suppress=True, precision=5)
        doctest.testmod()
        return 0
    if settings.test:
        settings.name = 'Kinetoplast'
        sequence = Sequence.KINETOPLAST
    elif sequence:
        settings.name = settings.name[:16]
        sequence = ''.join(sequence)
    else:
        parser.error('no sequence specified')

    def file_ext(f, e):
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
            results.save_csv(settings.csv)
        if settings.pdb:
            results.save_pdb(settings.pdb)
        if settings.seq:
            results.sequence.save(settings.seq)
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


MODELS = sorted(
    (a for a in dir(Model) if not a.startswith('_') and a.isupper()),
    key=lambda x: getattr(Model, x)['name']
)

if __name__ == '__main__':
    sys.exit(main())
