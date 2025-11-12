"""
Driver code for processing files and generating plots.

MIT License
-----------

Copyright (c) 2025 Daryl Gohl, Paul Egeler

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import subprocess
from enum import Enum
from types import SimpleNamespace

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import pysam


class OutputFormat(Enum):
    eps = 'eps'
    jpeg = 'jpeg'
    jpg = 'jpg'
    pdf = 'pdf'
    pgf = 'pgf'
    png = 'png'
    ps = 'ps'
    raw = 'raw'
    rgba = 'rgba'
    svg = 'svg'
    svgz = 'svgz'
    tif = 'tif'
    tiff = 'tiff'


def concatmap(args: SimpleNamespace) -> None:
    # TODO: Change signature
    ...
