"""
Algorithm by Robert Bridson.
"""


from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
import itertools

pi = np.pi


class PoissonDiskContext():
    """
    """

    def __init__(self, domain: np.ndarray, r) -> PoissonDiskContext:
        self.domain = domain
        self.n = len(domain)
        self.r = r
        self.cell_size = np.full(self.n, self.r / np.sqrt(self.n))
        self.grid_dims = np.array([np.ceil(length / cell_dim)
                                   for length, cell_dim in zip(self.domain, self.cell_size)]).astype(int)
        self.index_grid = np.full(self.grid_dims, -1)

        self.active_list = []
        self.point_list = []
        self.current_index = 0

    def __repr__(self):
        pass

    def which_cell(self, point) -> np.ndarray:
        """
        """

        point_cell = [int(np.floor(coord / grid_cell_dim))
                      for coord, grid_cell_dim in zip(point, self.cell_size)]
        return np.array(point_cell)

    def sample_from_annulus(self, point: np.ndarray, r1, r2) -> np.ndarray:
        """
        """

        n = len(point)

        # https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates

        r = np.random.uniform(r1, r2)
        # n-1 angles required to define a point in n dimensional space
        angles = [np.random.uniform(0, pi) for _ in range(n - 2)]
        angles.append(np.random.uniform(0, 2 * pi))

        if n == 2:
            offset_x = r * np.cos(angles[0])
            offset_y = r * np.sin(angles[0])
            offset = np.array([offset_x, offset_y])
        elif n == 3:
            offset_x = r * np.sin(angles[0]) * np.cos(angles[1])
            offset_y = r * np.sin(angles[0]) * np.sin(angles[1])
            offset_z = r * np.cos(angles[0])
            offset = np.array([offset_x, offset_y, offset_z])
        else:
            offset = np.array([r * np.sin(angles[:i]).prod() *
                               np.cos(angles[i]) for i in range(len(angles))])

        print(f"angles={angles}")
        print(f"offset={offset}")
        new_coords = point + offset

        return new_coords

    def conflicts(self, point1, point2, d) -> bool:
        """
        """

        assert point1.shape == point2.shape
        n = len(point1)

        dist = np.sqrt(((point1 - point2)**2).sum())

        if dist > d:
            return False
        else:
            return True

    def find_neighbors(self, cell):
        """
        """

        n = len(cell)
        offset_tuples = list(itertools.product(range(-1, 2), repeat=n))

        neighboring_indices = [tuple(cell + np.array(tup))
                               for tup in offset_tuples]

        return neighboring_indices

    def pd_sample(self, max_tries):
        """
        """

        k = 0
        while k < max_tries:

            # pick a point from the active list
            sel_from_active = np.random.randint(0, len(self.active_list))
            i = self.active_list[sel_from_active]
            start_point = self.point_list[i]

            # find a new point within the annulus of the selected point
            new_point = self.sample_from_annulus(
                start_point, self.r, 2 * self.r)
            test = [(point_coord > dom_lim) or (point_coord < 0)
                    for point_coord, dom_lim in zip(new_point, self.domain)]
            while np.any(test):
                new_point = self.sample_from_annulus(
                    start_point, self.r, 2 * self.r)  # resample
                test = [(point_coord > dom_lim) or (point_coord < 0)
                        for point_coord, dom_lim in zip(new_point, self.domain)]

            # find which cell the new point is in
            new_point_cell = self.which_cell(new_point)

            # find indices of neighboring cells
            neighboring_indices = self.find_neighbors(new_point_cell)

            # check non-empty neighboring cells for conflicts with new_point
            conflict_list = []
            for neighbor in neighboring_indices:
                try:
                    if self.index_grid[neighbor] != -1:
                        neighbor_point = self.point_list[self.index_grid[neighbor]]
                        print("I checked for a conflict")
                        c = self.conflicts(neighbor_point, new_point, self.r)
                        conflict_list.append(c)
                except:
                    conflict_list.append(False)

            if np.array(conflict_list).sum() == 0:
                # accept the new point
                print("I accept")
                return tuple(new_point_cell), new_point

            k += 1

        # return None if no suitable points were found
        return None, sel_from_active

    def _make_first_point(self) -> None:
        """
        """

        first_point = np.array([np.random.uniform(0, dim)
                                for dim in self.domain])
        point_cell = self.which_cell(first_point)
        self.index_grid[tuple(point_cell)] = self.current_index

        self.active_list.append(self.current_index)
        self.point_list.append(first_point)

    def generate(self) -> PoissonDiskContext:
        """
        """

        self._make_first_point()

        while self.active_list != []:
            c, s = self.pd_sample(30)

            if c is not None:
                self.current_index += 1
                self.point_list.append(s)
                self.index_grid[c] = self.current_index
                self.active_list.append(self.current_index)
            else:
                del self.active_list[s]

        return self
