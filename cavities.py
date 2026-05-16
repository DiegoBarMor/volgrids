import numpy as np

def display(array, char: str = 'x'):
    def render(c) -> str:
        if isinstance(c, (bool, np.bool)): return char if c else ' '
        return str(c)[0]

    print("-----------", array.shape)
    print('\n'.join(
        ''.join(render(c) for c in row)
        for row in array
    ))


tiles = (
    "                        ",
    "  xx                    ",
    " xxxx  xxxxxxxxxx xxx   ",
    "  xx  xxxxxx    xxxxxx  ",
    "  xx xxxxxx       xxx   ",
    "   xxxxxxxxxxxx         ",
    "  xxxx   xxxxxxxx       ",
    "   xxxx   xxxxxx        ",
    " xxxxxxxxxxxxxx         ",
    "    xxxxxxxxxx          ",
    "        xxxx            ",
    "                        ",
)

arr = np.array([
    [c == 'x' for c in row] for row in tiles
], dtype = bool)
display(arr)


vsurf = np.zeros_like(arr, dtype = bool)
hsurf = np.zeros_like(arr, dtype = bool)
vsurf[1:,:] = arr[1:,:] ^ arr[:-1,:]
hsurf[:,1:] = arr[:,1:] ^ arr[:,:-1]
surface = vsurf | hsurf
# display(surface)
# display(vsurf, '-')
# display(hsurf, '|')

mat_vrange = np.broadcast_to(np.arange(arr.shape[0])[:,None], arr.shape)
mat_hrange = np.broadcast_to(np.arange(arr.shape[1])[None,:], arr.shape)
# display(mat_vrange)

# vsurf_min = np.argmin(vsurf, axis = 0)
# hsurf_min = np.argmin(hsurf, axis = 1)
# print(vsurf_min)


vsurf_start = np.argmax(vsurf, axis = 0)
hsurf_start = np.argmax(hsurf, axis = 1)
vsurf_end   = np.argmax(vsurf[::-1,:], axis = 0)
hsurf_end   = np.argmax(hsurf[:,::-1], axis = 1)
# print(vsurf_start)
print(vsurf_end)
print(mat_vrange.shape[0] - vsurf_end - 1)

print(mat_vrange)

# mat_available = np.where(surface)
# display(mat_available)
vsurf_start[vsurf_start == 0] = arr.shape[0] # 0 implies no surface, so set the volume start to end of the dimension
hsurf_start[hsurf_start == 0] = arr.shape[1]

mat_available_top = mat_vrange < vsurf_start
# mat_available_bottom = mat_vrange >  vsurf_end
mat_available_bottom = mat_vrange >= (mat_vrange.shape[0] - vsurf_end - 1)
mat_available_left = (mat_hrange.T < hsurf_start).T
mat_available_right = (mat_hrange.T > (mat_hrange.shape[1] - hsurf_end)).T
# display(mat_available_top)
# display(mat_available_left)
display(mat_available_bottom)
# display(mat_available_right)
# # mat_available = np.less(mat_vrange, vsurf_start)
# # display(mat_available_top)
# # display(mat_available_bottom)
# display(mat_available_top | mat_available_bottom, char = '.')
# display(mat_available_left | mat_available_right, char = '.')
mat_av = mat_available_top.astype(int) + mat_available_bottom.astype(int) + mat_available_left.astype(int) + mat_available_right.astype(int)
cavities = 2 - mat_av
cavities[arr] = 0

# # display(mat_av)
# # display(surface)
# display(cavities)
