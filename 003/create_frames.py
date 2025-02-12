import pynbody
import matplotlib.pyplot as plt

IMAGE_PATH = "/home/lucasbondep/ic_astronomia/main/003/frames/"
SNAPSHOT_PATH = "/mnt/d/IC/snapshots/"

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

snapshots = []
lines = open(SNAPSHOT_PATH + "snapshot.txt", "r").readlines()
for i in range(len(lines)):
    file = f"{SNAPSHOT_PATH}snapshot_{i:03d}.hdf5"
    snapshots.append(file)

count = 0
for file in snapshots:
    print(f"{bcolors.OKGREEN}Generating frame {count} of {len(lines) - 1:03d}...{bcolors.ENDC}")
    
    s = pynbody.load(file)
    s.physical_units()

    mi = 0.6 # Average molecular weight
    Mh = 1.67262192 * 10**(-27) # Proton mass in kg
    s.gas["kT"] = (s.gas["u"] * (2 * mi * Mh) / 3) * 6.241506 * 10**15 * 10**(6)

    pynbody.plot.image(s.dm, qty="rho", width=8000, cmap="twilight", vmin=1e-9, vmax=1e9)
    plt.savefig(f"{IMAGE_PATH}dm/dm_{count:03d}.png", bbox_inches="tight", dpi=300)

    pynbody.plot.image(s.gas, qty="rho", width=8000, cmap="inferno", vmin=1e0, vmax=1e7)
    plt.savefig(f"{IMAGE_PATH}gas/gas_{count:03d}.png", bbox_inches="tight", dpi=300)

    pynbody.plot.image(s.gas, qty="kT", width=6000, cmap="inferno", log=False, vmin=0, vmax=20)
    plt.savefig(f"{IMAGE_PATH}gas/temp_{count:03d}.png", bbox_inches="tight", dpi=300)

    count += 1
