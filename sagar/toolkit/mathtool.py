# -*- coding: utf-8 -*-
import numpy

from math import sqrt
from itertools import combinations

def is_int_np_array(npa, atol=1e-5):
    return numpy.all(numpy.isclose(npa, numpy.around(npa), atol=atol))


def refine_positions(npa, atol=1e-5):
    """
    给定一个精度，让靠近1的元素变成零。则使得坐标在0～1之间可以直接比较。
    """
    flatten_npa = npa.flatten()
    for i, val in enumerate(flatten_npa):
        if val < atol or val > 1 - atol:
            flatten_npa[i] = 0
        else:
            flatten_npa[i] = val
    return flatten_npa.reshape(npa.shape)


def extended_gcd(aa, bb):
    """
    Algorithm: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Iterative_method_2

    parameters:
    aa, bb: int

    return: r, s, t

    r = s * aa + t * bb
    """
    lastremainder, remainder = abs(aa), abs(bb)
    x, lastx, y, lasty = 0, 1, 1, 0
    while remainder:
        lastremainder, (quotient, remainder) = remainder, divmod(
            lastremainder, remainder)
        x, lastx = lastx - quotient * x, x
        y, lasty = lasty - quotient * y, y
    return lastremainder, lastx * (-1 if aa < 0 else 1), lasty * (-1 if bb < 0 else 1)

def binomialCoeff(n, k):
    if k < 0:
        return -1
    result = 1
    if n - k < k:
        k = n - k

    for i in range(1, k + 1):
        result = result * (n - i + 1) / i
    return int(result)


def distance(p1, p2):
    #return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)
    p1, p2 = numpy.array(p1), numpy.array(p2)
    return numpy.linalg.norm(p1 - p2)

def closest_pair(points):
    """
    points is a list or numpy.array of 1x3 cartesian positions
    algorithm: http://en.wikipedia.org/wiki/Closest_pair_of_points
    (brute force) O(n^2)
    """
    # TODO: using O(n log n) algorithm
    num_points = len(points)
    if num_points < 2:
        return float('inf')
    return min([distance(points[i], points[j])
                for i, j in combinations(range(num_points), 2)])
