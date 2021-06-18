#!/usr/bin/env python3
# dnacurve_web.py

# Copyright (c) 2005-2021, Christoph Gohlke
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

"""DNA Curvature Analysis - A simple Web interface for dnacurve.py.

Run ``python -m dnacurve --web`` to execute the script in a local web server.

Revisions
---------
2021.6.18
    Link to source code on GitHub.
    Fix failures on WSL2.
2021.3.6
    Update copyright and formatting.
2020.1.1
    Fix csv file extension.
    Remove support for Python 2.7 and 3.5.
    Update copyright.
2018.8.15
    Move module into dnacurve package.
2018.5.29
    Use CSS flex layout.
    Separate styles from content.
    Run http server from this file's directory.
2018.5.25
    Backwards incompatible changes:
    Rename to dnacurve_web.py.
    Accept Flask request.form.
    Style and template changes.
    Display models inline.
2018.3.11
    Use Data URLs in HTML instead of external files.
2018.2.6
    Style fixes.
    Escape error messages.
2005.x.x
    Initial release.

"""

__version__ = '2021.6.18'

import os
import sys
import io
import base64
import hashlib
from html import escape

try:
    from . import dnacurve
except ImportError:
    import dnacurve


PAGE = """<!DOCTYPE html PUBLIC
"-//W3C//DTD XHTML 1.1 plus MathML 2.0 plus SVG 1.1//EN"
"http://www.w3.org/2002/04/xhtml-math-svg/xhtml-math-svg.dtd">
<html xmlns="http://www.w3.org/1999/xhtml"
xmlns:mathml="http://www.w3.org/1998/Math/MathML"
xmlns:svg="http://www.w3.org/2000/svg" xml:lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<style type="text/css">
html {{font-size:16px}}
body {{min-width:576px}}
div.header {{display:flex;flex-wrap:wrap;align-items:baseline}}
div.header a {{font-weight:bold;text-decoration:none}}
h1 {{margin:0; padding:0 0.5em 0 0}}
h1 a {{text-decoration:none;color:#000000}}
form {{display:flex;flex-direction:column;flex-wrap:wrap;
background-color:#eeeeee;border:1px solid #aaaaaa;padding:1em 1em 0 1em}}
form a {{font-weight:bold;margin-left:0.3em}}
input {{margin-left:0.5em;min-width:8em}}
label {{font-weight:bold}}
textarea {{width:100%;height:8em}}
div.row {{display:flex;flex-wrap:wrap;justify-content:space-between;
align-items:center;margin-top:1em}}
div.row div {{margin-bottom:1em}}
div.buttons {{align-self:flex-end;margin-bottom:1em;margin-left:auto}}
div.content img {{width:576px;height:720px;padding:0}}
div.content pre {{font-size:75%}}
</style>
{heads}
<meta name="generator" content="dnacurve_web.py" />
<meta name="robots" content="noarchive" />
<meta name="format-detection" content="telephone=no" />
<meta name="viewport" content="width=608px" />
<title>DNA Curvature Analysis - Christoph Gohlke</title>
</head>
<body>
<div class="header">
<h1><a href="{url}"
title="DNA Curvature Analysis version {version} by Christoph Gohlke"
>DNA Curvature Analysis</a></h1>
<p>by <a href="https://www.lfd.uci.edu/~gohlke/">Christoph Gohlke</a></p>
</div>
<form id="dnacurve" method="post" action="">
<div>
<label for="seq">Sequence:</label>
<span>(max {maxlen:d} nucleotides)</span>
</div>
<div class="textarea">
<textarea name="seq" id="seq">{sequence}</textarea>
</div>
<div class="row">
<div>
<label for="mod">Model:</label>
<select name="mod" id="mod">
{models}
</select>
<a href="javascript:document.forms.dnacurve.q.value='models';
document.forms.dnacurve.submit()" rel="nofollow" title="Show models">?</a>
</div>
<div class="buttons">
<input type="reset" value="Reset" onclick="window.location='{url}'"/>
<input type="submit" value="Submit" />
<input type="hidden" name="q" />
</div>
</div>
</form>
<div class="content">
{content}
</div>
</body>
</html>"""

HELP = """<p>This web application calculates the global 3D structure of a
DNA molecule from its nucleotide sequence according to the dinucleotide
wedge model.
Local bending angles and macroscopic curvature are analyzed.</p>
<p>Try the
<a href="" onclick="javascript:document.forms.dnacurve.seq.value='{s1}';
return false;">Kinetoplast</a> or
<a href="" onclick="javascript:document.forms.dnacurve.seq.value='{s2}';
return false;">Phased AAAAAA</a> sequences for example.</p>
<p>For each nucleotide at position i of the input sequence,
the following values are calculated:</p>
<ul>
<li>the <strong>3D coordinates</strong> of the helix axis and
5' Phosphate atoms of a B-DNA.</li>
<li>the <strong>curvature</strong>, which is the inverse of
the radius of a circle passing through helix axis coordinates at
i-10, i, and i+10, relative to the curvature in a nucleosome.</li>
<li>the <strong>curvature angle</strong> between the smoothed
basepair normal vectors at i-15 and i+15, relative to the curvature
in a nucleosome.</li>
<li>the local <strong>bend angle</strong> between the basepair
normal vectors at i-2 and i+2.</li>
</ul>
<p><strong>Output files</strong> of the DNA curvature analysis are:</p>
<ul>
<li>a plot of the helix backbone coordinates and curvature values.</li>
<li>all calculated values in CSV format.</li>
<li>the helix backbone coordinates in PDB format.</li>
</ul>
<div class="references">
<h3>References</h3>
<ul>
<li><a href="https://www.ncbi.nlm.nih.gov/pubmed/2006170">Curved DNA
without A-A: experimental estimation of all 16 DNA wedge angles</a>.
<br />Bolshoy A. et al; Proc Natl Acad Sci USA 88; 2312-6 (1991)</li>
<li><a href="https://www.ncbi.nlm.nih.gov/pubmed/7816643">Bending and
curvature calculations in B-DNA</a>.<br />Goodsell DS; Dickerson RE;
Nucleic Acids Res 22; 5497-503 (1994)</li>
<li><a href="https://www.ncbi.nlm.nih.gov/pubmed/3271483">A comparison
of six DNA bending models</a>.<br />Tan RK; Harvey SC;
J Biomol Struct Dyn 5; 497-512 (1987)</li>
<li><a href="https://www.ncbi.nlm.nih.gov/pubmed/3456570">Curved DNA:
design, synthesis, and circularization</a>.<br />Ulanovsky L. et al;
Proc Natl Acad Sci USA 83; 862-6 (1986)</li>
</ul>
</div>
<div class="disclaimer">
<h3>Disclaimer</h3>
<p>Because this service is provided free of charge, there is no
warranty for the service, to the extent permitted by applicable law.
The service is provided &quot;as is&quot; without warranty of any kind,
either expressed or implied, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose.
The entire risk as to the quality and performance is with you.</p>
</div>
<div class="about">
<h3>About</h3>
<p>DNA Curvature Analysis version {version} by
<a href="https://www.lfd.uci.edu/~gohlke/">Christoph Gohlke</a>.
Source code is available on
<a href="https://github.com/cgohlke/dnacurve">GitHub</a>.</p>
</div>"""

RESULTS = """<!--h2>Results for {fname}</h2-->
<img src="data:{mime};name={fname}.{imageext};base64,{plot}" alt="Plot" />
<ul>
<li><a href='data:text/csv;name={fname}.csv;base64,{csv}'
download='{fname}.csv'>Calculated values</a> (CSV format)</li>
<li><a href='data:chemical/x-pdb;name={fname}.pdb;base64,{pdb}'
download='{fname}.pdb'>Helix coordinates</a> (PDB format)</li>
</ul>"""


def response(
    form,
    url,
    template=PAGE,
    help=HELP,
    maxlen=dnacurve.MAXLEN,
    analyze_function=None,
    heads='',
):
    """Return HTML document from submitted form data."""
    seq = form.get('seq', '').strip()
    mod = form.get('mod', '')

    if form.get('q', '') == 'version':
        content = f'<p>Version: {dnacurve.__version__}</p>'
    elif form.get('q', '') == 'models':
        content = ['<h2>Curvature Models</h2>']
        for model in dnacurve.MODELS:
            lines = str(dnacurve.Model(model)).splitlines()
            content.append(f'<h3>{escape(lines[0])}</h3>')
            content.append(
                '<pre>{}</pre>'.format(escape('\n'.join(lines[1:])))
            )
        content = '\n'.join(content)
    elif seq:
        if analyze_function is None:
            analyze_function = analyze
        content = analyze_function(seq, mod, maxlen)
    else:
        content = help.format(
            version=dnacurve.__version__,
            s1=''.join(dnacurve.Sequence.KINETOPLAST.split())[:maxlen],
            s2=(dnacurve.Sequence.PHASED_AAAAAA * 14)[:maxlen],
        )

    options = []
    for model in dnacurve.MODELS:
        if model == mod:
            option = '<option value="{}" selected="selected">{}</option>'
        else:
            option = '<option value="{}">{}</option>'
        label = getattr(dnacurve.Model, model)['name']
        options.append(option.format(escape(model), escape(label)))
    options = '\n'.join(options)

    return template.format(
        sequence=escape(seq),
        models=options,
        content=content,
        url=url,
        version=dnacurve.__version__,
        heads=heads.strip(),
        maxlen=maxlen,
    )


def analyze(sequence, model, maxlen, template=RESULTS, imageformat='svg'):
    """Return results of DNA curvature analysis as HTML string."""
    try:
        seq = dnacurve.Sequence(sequence)
        if len(seq) > maxlen:
            raise ValueError('sequence is too long')
        if len(seq) < 10:
            raise ValueError('sequence is too short')
        name = hashlib.md5(seq.string.encode('ascii')).hexdigest()
        seq.name = name[:13]
        fname = name + str(dnacurve.MODELS.index(model))
        mod = dnacurve.Model(model)
        results = dnacurve.CurvedDNA(seq, mod)
        plot = io.BytesIO()
        results.plot(plot, imageformat=imageformat)
        plot = b64encode(plot.getvalue())
        csv = b64encode(results.csv())
        pdb = b64encode(results.pdb())
    except Exception as e:
        # raise
        e = str(e).splitlines()
        text = e[0][0].upper() + e[0][1:]
        details = '\n'.join(e[1:])
        return '<h2>Error</h2><p>{}</p><pre>{}</pre>'.format(
            escape(text), escape(details)
        )
    mime = {'svg': 'image/svg+xml', 'png': 'image/png'}[imageformat]
    return template.format(
        plot=plot,
        pdb=pdb,
        csv=csv,
        fname=fname,
        mime=mime,
        imageext=imageformat,
    )


def b64encode(s):
    """Return argument as b64encode encoded str."""
    try:
        s = s.encode('ascii')
    except Exception:
        pass
    return base64.b64encode(s).decode(errors='strict')


def main(url='http://localhost:9000/{}'):
    """Run web application in local web server."""
    import cgi
    import cgitb

    cgitb.enable()

    dirname, filename = os.path.split(__file__)
    if filename[-1] != 'y':
        filename = filename[:-1]  # don't use .pyc or .pyo
    url = url.format(filename)
    if dirname:
        os.chdir(dirname)

    if os.getenv('SERVER_NAME'):
        print('Content-type: text/html\n\n')
        request = cgi.FieldStorage()
        request.get = request.getfirst
        print(response(request, url))
    else:
        import webbrowser
        from urllib.parse import urlparse
        from http.server import HTTPServer, CGIHTTPRequestHandler

        def is_cgi(self):
            """Monkey patch for CGIHTTPRequestHandler.is_cgi."""
            if filename in self.path:
                self.cgi_info = '', self.path[1:]
                return True
            return False

        CGIHTTPRequestHandler.is_cgi = is_cgi
        print('Serving CGI script at', url)
        webbrowser.open(url)
        url = urlparse(url)
        HTTPServer(
            (url.hostname, url.port), CGIHTTPRequestHandler
        ).serve_forever()


if __name__ == '__main__':
    sys.exit(main())
