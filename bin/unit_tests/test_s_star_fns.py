import s_star_fns
import numpy as np
from numpy.testing import assert_almost_equal as aae
import sys


def test_calc_geno_dist():
    def old_calc_geno_dist(gt1, gt2):
        gd = [abs(gt1[i] - gt2[i]) for i in range(len(gt1))]
        return sum(gd)
    np.random.seed(0)
    for _ in range(20):
        gt1 = np.random.rand(10)
        gt2 = np.random.rand(10)
        gt1list = gt1.tolist()
        gt2list = gt2.tolist()

        # old version requires list, new takes either
        aae(s_star_fns.calc_geno_dist(gt1, gt2),
            old_calc_geno_dist(gt1list, gt2list))
        aae(s_star_fns.calc_geno_dist(gt1list, gt2list),
            old_calc_geno_dist(gt1list, gt2list))

    assert s_star_fns.calc_geno_dist([], []) == 0


def test_calc_s():
    bonus = 10
    mismatch = 3
    penalty = -5

    assert s_star_fns.calc_s([], [],
                             1, 2,
                             bonus, mismatch, penalty) == -sys.maxint

    assert s_star_fns.calc_s([], [],
                             11, 2,
                             bonus, mismatch, penalty) == -sys.maxint

    assert s_star_fns.calc_s([], [],  # gd == 0
                             12, 2,
                             bonus, mismatch, penalty) == bonus + 10  # dist

    assert s_star_fns.calc_s([1], [0],  # gd == 1
                             12, 2,
                             bonus, mismatch, penalty) == penalty

    assert s_star_fns.calc_s([1, 1, 1], [0, 0, 0],  # gd == 3
                             12, 2,
                             bonus, mismatch, penalty) == penalty

    assert s_star_fns.calc_s([1, 1, 1, 0], [0, 0, 0, 1],  # gd == 4
                             12, 2,
                             bonus, mismatch, penalty) == -sys.maxint


def test_calc_s_star():
    print(s_star_fns.calc_s_star(
        [1, 1, 2, 1],
        [10, 15, 30, 45],
        5000, 5, -10000))
    print(s_star_fns.calc_s_star(
        [1, 1, 1, 1],
        [5382, 28610, 32662, 35985],
        5000, 5, -10000))
    # TODO finish this assert 0
