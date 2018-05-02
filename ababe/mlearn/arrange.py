# coding: utf-8
# Distributed under the terms of the MIT License.
import ababe.mlearn.dsutils as dsutils

class ArrangeProblem(object):

    def get_start_position(self):
        dsutils.raise_not_defined()

    def is_goal(self, tmp_result):
        dsutils.raise_not_defined()

    def get_neighbors(self, state):
        dsutils.raise_not_defined()

def random_walk_arrange(problem):
    path = []
    start = problem.get_start_position()
    path.append(start)
    # frontier is a random pop DataStur should be a set()
    frontier = dsutils.RandomPopSet()

    for pos in problem.get_neighbors(start):
        frontier.push(pos)

    while not problem.is_goal(len(path)):
        current = frontier.random_pop()
        path.append(current)

        for pos in problem.get_neighbors(current):
            if pos not in path:
                frontier.push(pos)

    # path is a list contain the sequence of position
    return path