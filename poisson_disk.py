import numpy as np
from matplotlib import pyplot as plt

pi = np.pi

def which_cell(point, grid_cell_size) -> np.ndarray:
    """
    """

    point_cell = [int(np.floor(coord / grid_cell_dim)) for coord, grid_cell_dim in zip(point, grid_cell_size)]
    return np.array(point_cell)

def sample_from_annulus(point : np.ndarray, r1, r2) -> np.ndarray:
    """
    """

    n = len(point)

    r = np.random.uniform(r1, r2)
    angles = [np.random.uniform(0, pi) for _ in range(n-2)] # n-1 angles required to define a point in n dimensional space
    angles.append(np.random.uniform(0, 2*pi))

    if n == 2:
        offset_x = r*np.cos(angles[0])
        offset_y = r*np.sin(angles[0])
        offset = np.array([offset_x, offset_y])
    else:
        offset = np.array([r * np.sin(angles[:i]).prod() * np.cos(angles[i]) for i in range(len(angles))])

    new_coords = point + offset

    return new_coords

def conflicts(point1, point2, d) -> bool:
    """
    """

    assert point1.shape == point2.shape
    n = len(point1)

    dist = np.sqrt(((point1 - point2)**2).sum())

    if dist > d:
        return False
    else:
        return True




def pd_sample(min_dist, max_tries, domain, active_list, point_list, cell_size, index_grid):
    """
    """

    k = 0
    while k < max_tries:

        print(f"k={k}")

        # pick a point from the active list
        sel_from_active = np.random.randint(0, len(active_list))
        i = active_list[sel_from_active]
        start_point = point_list[i]

        # find a new point within the annulus of the selected point
        new_point = sample_from_annulus(start_point, min_dist, 2*min_dist)
        test = [(point_coord > dom_lim) or (point_coord < 0) for point_coord, dom_lim in zip(new_point, domain)]
        while np.any(test):
            new_point = sample_from_annulus(start_point, min_dist, 2*min_dist) # resample
            test = [(point_coord > dom_lim) or (point_coord < 0) for point_coord, dom_lim in zip(new_point, domain)]

        # find which cell the new point is in
        new_point_cell = which_cell(new_point, cell_size)

        # find indices of neighboring cells
        neighboring_indices = [new_point_cell + np.array([x, y]) for x in range(-1, 2) for y in range(-1, 2)]
        neighboring_indices = [tuple(index) for index in neighboring_indices]

        # check non-empty neighboring cells for conflicts with new_point
        conflict_list = []
        for neighbor in neighboring_indices:
            try:
                if index_grid[neighbor] != -1:
                    neighbor_point = point_list[index_grid[neighbor]]
                    print("I checked for a conflict")
                    c = conflicts(neighbor_point, new_point, min_dist)
                    conflict_list.append(c)
            except: # this is ALWAYS going to the except block for some reason
               conflict_list.append(False)


        if np.array(conflict_list).sum() == 0:
            # accept the new point
            print("I accept")
            return tuple(new_point_cell), new_point

        k += 1

    # return None if no suitable points were found
    return None, sel_from_active
