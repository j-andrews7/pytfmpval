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

See full documentation and use examples `here <http://pytfmpval.readthedocs.io/en/latest/>_`.

Installation
---------------

pytfmpval is on PyPI, so you can install via ``pip`` easily::

    pip install pytfmpval


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