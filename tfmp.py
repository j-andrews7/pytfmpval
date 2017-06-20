#!/usr/bin/env python
"""
Testing pytfmpval functions.
"""

import pytfmpval as tfm
from math import ceil


def create_matrix(matrix_file, bg=[0.25, 0.25, 0.25, 0.25], mat_type="pfm", log_type="nat"):
    """
    From a JASPAR formatted matrix file, create a Matrix object.

    This function also converts it to a log-odds (position weight) matrix if necessary.

    args:
        matrix_file (str):
            Name of file containing matrix.
        bg (list of floats):
            Background nucleotide frequencies for [A, C, G, T].

    returns:
        m (pytfmpval Matrix):
            Matrix in pwm format.
    """

    a, c, g, t = bg[0], bg[1], bg[2], bg[3]

    m = tfm.Matrix(a, c, g, t)
    m.readJasparMatrix(matrix_file)

    if mat_type == "pfm":
        if log_type == "nat":
            m.toLogOddRatio()
        else:
            m.toLog2OddRatio()

    return m


def read_matrix(matrix, bg=[0.25, 0.25, 0.25, 0.25], mat_type="pfm", log_type="nat"):
    """
    From a JASPAR formatted matrix file, create a Matrix object.

    This function also converts it to a log-odds (position weight) matrix if necessary.

    args:
        matrix_file (str):
            Name of file containing matrix.
        bg (list of floats):
            Background nucleotide frequencies for [A, C, G, T].

    returns:
        m (pytfmpval Matrix):
            Matrix in pwm format.
    """

    a, c, g, t = bg[0], bg[1], bg[2], bg[3]

    m = tfm.Matrix(a, c, g, t)
    m.readMatrix(matrix)

    if mat_type == "pfm":
        if log_type == "nat":
            m.toLogOddRatio()
        else:
            m.toLog2OddRatio()

    return m


def score2pval(matrix, req_score):
    """
    Determine the p-value for a given score for a specific motif PWM.

    args:
        matrix (pytfmpval Matrix):
            Matrix in pwm format.
        req_score (float):
            Requested score for which to determine the p-value.

    returns:
        ppv (float):
            The calculated p-value corresponding to the score.
    """

    granularity = 0.1
    max_granularity = 1e-10
    decrgr = 10  # Factor to increase granularity by after each iteration.

    pv = tfm.doublep()
    ppv = tfm.doublep()

    while granularity > max_granularity:
        matrix.computesIntegerMatrix(granularity)
        max_s = int(req_score * matrix.granularity + matrix.offset + matrix.errorMax + 1)
        min_s = int(req_score * matrix.granularity + matrix.offset - matrix.errorMax - 1)
        score = int(req_score * matrix.granularity + matrix.offset)

        matrix.lookForPvalue(score, min_s, max_s, ppv, pv)

        if ppv.value() == pv.value():
            return ppv.value()

        granularity = granularity / decrgr

    print("Max granularity exceeded. Returning closest approximation.")
    return ppv.value()


def pval2score(matrix, pval):
    """
    Determine the score for a given p-value for a specific motif PWM.

    args:
        matrix (pytfmpval Matrix):
            Matrix in pwm format.
        pval (float):
            p-value for which to determine the score.

    returns:
        score (float):
            The calculated score corresponding to the p-value.
    """

    init_granularity = 0.1
    max_granularity = 1e-10
    decrgr = 10  # Factor to increase granularity by after each iteration.

    pv = tfm.doublep()  # Initialize as a c++ double.
    ppv = tfm.doublep()
    matrix.computesIntegerMatrix(init_granularity)
    max_s = matrix.maxScore + ceil(matrix.errorMax + 0.5)
    min_s = matrix.minScore
    granularity = init_granularity

    while granularity > max_granularity:
        matrix.computesIntegerMatrix(granularity)

        score = matrix.lookForScore(min_s, max_s, pval, pv, ppv)

        min_s = (score - ceil(matrix.errorMax + 0.5)) * decrgr
        max_s = (score + ceil(matrix.errorMax + 0.5)) * decrgr

        if ppv.value() == pv.value():
            break

        granularity = granularity / decrgr

    if granularity <= max_granularity:
        print("Max granularity exceeded. Returning closest score approximation.")

    final_score = (score - matrix.offset) / matrix.granularity
    return final_score
