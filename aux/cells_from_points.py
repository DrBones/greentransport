def cells_from_points(array):
    from scipy.ndimage.interpolation import geometric_transform
    def shift_func(output_coords):
        return (output_coords[0] - 0.5, output_coords[1] - 0.5)
    cells = geometric_transform(array, shift_func, output_shape=(array.shape[0]-1,array.shape[1]-1))
    return cells
