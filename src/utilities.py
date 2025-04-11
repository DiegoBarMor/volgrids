import json, time
import numpy as np
from scipy import interpolate

################################################################################ CLASSES
class Timer:
    def __init__(self, display = ''):
        if display: print(display, end = ' ', flush = True)
        self.start = time.time()

    def end(self):
        elapsed = time.time() - self.start
        print(f"({int(elapsed // 60)}m {elapsed % 60:.2f}s)", flush = True)


################################################################################ FUNCTIONS
def read_json(path_json):
    with open(path_json, 'r') as file:
        return json.load(file)

def write_json(path_json, content):
    with open(path_json, 'w') as file:
        file.write(json.dumps(content))

def save_metadata(metadata):
    meta = metadata.copy()
    for k in ("pdb", "out", "apbs", "meta"):
        meta[k] = str(meta[k])
    write_json(meta["meta"], meta)

def normalize(vector):
    """
    input:  (3,)
    output: (3,)
    """
    return vector / np.linalg.norm(vector)

def get_norm(m_vectors):
    """
    input:  (xres, yres, zres, 3)
    output: (xres, yres, zres)
    """
    return np.linalg.norm(m_vectors, axis = -1).T

def dot_product(m_vectors, vector):
    """
    input (m_vectors): (xres, yres, zres, 3)
    input (vector):    (3,)
    output:            (xres, yres, zres)
    """
    d = m_vectors[:,:,:,0] * vector[0] + \
        m_vectors[:,:,:,1] * vector[1] + \
        m_vectors[:,:,:,2] * vector[2]
    return d

def get_angle(m_vectors, vector, isStacking):
    """
    input (m): (xres, yres, zres, 3)
    input (v): (3,)
    output   : (xres, yres, zres)
    """
    numerator = dot_product(m_vectors, vector)
    denominator = get_norm(m_vectors).T * np.linalg.norm(vector).T

    mask = denominator == 0
    numerator[mask] = 1
    denominator[mask] = 1

    cos_val = np.clip(numerator / denominator, -1, 1)
    angle_radians = np.arccos(cos_val)
    angle_degrees = angle_radians * 180 / np.pi

    if isStacking:
        angle_degrees[angle_degrees >= 90] = 180 - angle_degrees[angle_degrees >= 90]
    else: # hydrogen bonds
        angle_degrees = 180 - angle_degrees

    angle_degrees[mask] = -90

    return angle_degrees.T

def univariate_gaussian(x, mu, sigma):
    """ input_mat "x" shape: (xsize, ysize, zsize), values: dist """
    u = x - mu
    s = 1 / (sigma ** 2)
    return np.exp(-(1/2) * np.power(u, 2) * s)

def multivariate_gaussian(x, mu, cov_inv):
    """ input_mat "x" shape: (xsize, ysize, zsize, 2 = (dist, beta)) """
    u = x - mu
    sigma = cov_inv

    ux, uy = u[:,:,:,0], u[:,:,:,1]
    a,b,c,d = sigma[0,0], sigma[0,1], sigma[1,0], sigma[1,1]
    return np.exp(-(1/2) * (ux * (a * ux + b * uy) + uy * (c * ux + d * uy)))

def interpolate_3d(x0, y0, z0, data_0, new_coords):
    return interpolate.RegularGridInterpolator(
        (x0, y0, z0), data_0, bounds_error = False, fill_value = 0
    )(new_coords).T

def format_vector_str(vector : np.array) -> str:
        return '(' + (
            ' '.join(vector.astype(str)) if vector.dtype == int else \
            ' '.join(map(lambda n: f"{n:.3f}", vector))
        ) + ')'


################################################################################
