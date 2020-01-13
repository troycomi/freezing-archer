import s_star_fns
import numpy as np
from numpy.testing import assert_almost_equal as aae
import sys


def old_calc_s(gt1, gt2, position1, position2,
               match_bonus=5000, max_mismatch=5, mismatch_penalty=-10000):

    # was abs, only use case with 1 < 2
    # dist = abs(position1 - position2)
    dist = position2 - position1
    if dist < 10:
        return np.NINF

    gd = old_calc_geno_dist([gt1], [gt2])
    if gd == 0:
        return match_bonus + dist

    elif gd <= max_mismatch:
        return mismatch_penalty

    return np.NINF


def old_calc_geno_dist(gt1, gt2):
    gd = [abs(gt1[i] - gt2[i]) for i in range(len(gt1))]
    return sum(gd)


def test_calc_s_dists():
    assert s_star_fns.calc_s_dists([], [], 5000, 5, -10000).size == 0

    bonus = 10
    mismatch = 3
    penalty = -5

    gen = [0, 1, 2, 0, 1, 2]
    pos = [1, 2, 15, 25, 45, 100]

    result = np.full((6, 6), 0.0)
    for i in range(6):
        for j in range(6):
            result[i, j] = old_calc_s(gen[i], gen[j], pos[i], pos[j],
                                      bonus, mismatch, penalty)
    aae(s_star_fns.calc_s_dists(gen, pos, bonus, mismatch, penalty),
        result)

    # more spots are ni
    mismatch = 1
    result = np.full((6, 6), 0.0)
    for i in range(6):
        for j in range(6):
            result[i, j] = old_calc_s(gen[i], gen[j], pos[i], pos[j],
                                      bonus, mismatch, penalty)
    aae(s_star_fns.calc_s_dists(gen, pos, bonus, mismatch, penalty),
        result)


def test_calc_s_star():
    print(s_star_fns.calc_s_star(
        [1, 1, 2, 1],
        [10, 15, 30, 45],
        5000, 5, -10000))
    print(s_star_fns.calc_s_star(
        [1, 1, 1, 1],
        [5382, 28610, 32662, 35985],
        5000, 5, -10000))
    # TODO assert 0
