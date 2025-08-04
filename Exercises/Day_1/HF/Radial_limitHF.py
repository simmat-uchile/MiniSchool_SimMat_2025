import numpy as np
import matplotlib.pyplot as plt

# === GRID PARAMETERS ===
delta_bohr = 0.03333 * 1.889726  # grid spacing in bohr
#delta_bohr = 0.03333  # grid spacing in bohr
voxel_volume = delta_bohr ** 3   # volume of one voxel
dr = 0.05                        # radial bin width (bohr)

# === FILES AND LABELS ===
files = [
    ("He_sto3g_density.xyz", "STO-3G"),
    ("He_augccpv6z_density.xyz", "aug-cc-pV6Z")
]

# === SMOOTHING FUNCTION ===
def moving_average(y, window=7):
    return np.convolve(y, np.ones(window)/window, mode='same')

# === PLOT SETUP ===
plt.figure()
rmax_global = 0.0  # to track the largest radius

# === PROCESS EACH FILE ===
for filename, label in files:
    data = np.loadtxt(filename)
    x, y, z, rho = data.T
    r = np.sqrt(x**2 + y**2 + z**2)

    rmax = np.max(r)
    rmax_global = max(rmax_global, rmax)

    nbins = int(np.ceil(rmax / dr))
    hist = np.zeros(nbins)
    counts = np.zeros(nbins)

    for i in range(len(r)):
        ir = int(r[i] / dr)
        if ir < nbins:
            hist[ir] += rho[i] * voxel_volume
            counts[ir] += 1

    rvals = (np.arange(nbins) + 0.5) * dr
    D_r = np.divide(hist, dr, out=np.zeros_like(hist), where=counts != 0)
    smooth_D_r = moving_average(D_r, window=7)

    total_electrons = np.sum(D_r) * dr
    print(f"{label} - Total number of electrons: {total_electrons:.6f}")

    plt.plot(rvals, smooth_D_r, label=label)

# === FINAL PLOT SETTINGS ===
plt.xlabel("r (bohr)")
plt.ylabel("D(r) (electrons / bohr)")
plt.title("Radial Distribution Function")
plt.legend()
plt.grid(False)
plt.tight_layout()
plt.xlim(0, rmax_global)  # set x-axis range based on data
plt.savefig("radial_distribution.png", dpi=300)

