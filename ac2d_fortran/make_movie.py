import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from scipy.io import FortranFile


amp = 1.
dt = 0.0005

NX = 301
NZ = 201

# Use if fortran output is single precision
# ftype = np.float32

# Use if fortran output is double precision
ftype = np.float64


def load_matrix(filename):
    f = FortranFile(filename, 'r')
    data = f.read_reals(ftype).reshape((NZ, NX), order="F")

    return data


def create_animation(files, show=True, savefn=None):
    data = load_matrix(files[0])
    print(data.shape)
    fig, ax = plt.subplots(1, 1)
    im = ax.imshow(data, cmap='seismic', vmin=-amp, vmax=amp)
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Z [m]")
    fig.colorbar(im)

    def update_image(i):
        data = load_matrix(files[i+1])
        im.set_data(data)
        ax.set_title(f"2D Acoustic Wave Propagation\nTime: {i*dt:10.5f} [s]")

    ani = FuncAnimation(fig, update_image, len(files)-1, interval=50)

    if show:
        plt.show()

    if savefn is not None:
        ani.save(savefn + ".mp4")

    return ani


def main():
    nx = 50
    nz = 50

    folder = "pressure"
    files = [os.path.join(folder, f) for f in sorted(os.listdir(folder))]

    create_animation(files)


if __name__ == "__main__":
    main()
