pytfmpval
===========

.. image:: https://travis-ci.org/j-andrews7/pytfmpval.png?branch=master
    :target: https://travis-ci.org/j-andrews7/pytfmpval

.. image:: https://badge.fury.io/py/pytfmpval.svg?style=flat
    :target: http://badge.fury.io/py/pytfmpval

.. image:: https://readthedocs.org/projects/pytfmpval/badge/?version=latest
    :target: http://pytfmpval.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

This Python package serves as a wrapper for the incredibly useful `TFM-Pvalue C++ program <http://bioinfo.lifl.fr/tfm-pvalue/tfm-pvalue.php>`_. It allows users to determine score thresholds for a given transcription factor position frequency matrix associated with a specific p-value. Naturally, it can also perform the reverse, quickly calculating an accurate p-value from a score for a given motif matrix.

``pytfmpval`` allows this functionality to be easily utilized within a Python script, module, or package.

See full documentation and use examples `here <http://pytfmpval.readthedocs.io/en/latest/>`_.

Installation
---------------

pytfmpval is on PyPI, so you can install via ``pip`` easily::

    pip install pytfmpval


A Simple Example
--------------------------

`JASPAR <http://jaspar.genereg.net>`_ is a very highly-touted transcription factor motif database from which motif count matrices can be downloaded for a large variety of organisms and transcription factors. There exist numerous other motif databases as well (TRANSFAC, CIS-BP, MEME, HOMER, WORMBASE, etc), most of which use a relatively similar format for their motifs. Typically, a motif file consists of four rows or columns with each position in a given row or column corresponding to a base within the motif. Sometimes there is an comment line started with ``>``. The row or column order is always ``A, C, G, T``. In this example, the motif consists of four rows corresponding to the 16 positions of the motif with counts for each base at each position.

>>> from pytfmpval import tfmp
>>> m = tfmp.create_matrix("MA0045.pfm")
>>> tfmp.score2pval(m, 8.7737)
9.992625564336777e-06
>>> tfmp.pval2score(m, 0.00001)
8.773708000000001

This could also be done by creating a string for the matrix by concatenating the rows (or columns) and using the ``read_matrix()`` function. This method is usually easier, as allows the user to parse the motif file as necessary to ensure a proper input. It's also more fitting for high-throughput.

>>> from pytfmpval import tfmp
>>> mat = (" 3  7  9  3 11 11 11  3  4  3  8  8  9  9 11  2" 
...        " 5  0  1  6  0  0  0  3  1  4  5  1  0  5  0  7"  
...    " 4  3  1  4  3  2  2  2  8  6  1  4  2  0  3  0" 
...    " 2  4  3  1  0  1  1  6  1  1  0  1  3  0  0  5"
...   )
>>> m = tfmp.read_matrix(mat)
>>> tfmp.pval2score(m, 0.00001)
8.773708000000001
>>> tfmp.score2pval(m, 8.7737)
9.992625564336777e-06

Contribute
---------------

Any and all contributions are welcome. Bug reporting via the `Issue Tracker <github.com/j-andrews7/pytfmpval/issues>`_ is much appeciated. Here's how to contribute:

1. Fork the `pytfmpval repository <https://github.com/j-andrews7/pytfmpval>`_ on github (see `forking help <https://help.github.com/articles/fork-a-repo/>`_).

2. Make your changes/fixes/improvements locally.

3. Optional, but much-appreciated: write some tests for your changes. (Don't worry about integrating your tests into the test framework - writing some in your commit comments or providing a test script is fine. I will integrate them later.)

4. Send a pull request (see `pull request help <https://help.github.com/articles/about-pull-requests/>`_).


Reference
--------------

| Efficient and accurate P-value computation for Position Weight Matrices
| H. Touzet and J.S. Varré
| *Algorithms for Molecular Biology 2007, 2:15*

License
-----------

This project is licensed under the GPL3 license. You are free to use, modify, and distribute it as you see fit. The program is provided as is, with no guarantees.