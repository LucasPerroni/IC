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
    print(f"{bcolors.OKGREEN}Generating frame {count}...{bcolors.ENDC}")
    
    s = pynbody.load(file)
    s.physical_units()

    pynbody.plot.image(s.dm, qty="rho", width=8000, cmap="twilight", vmin=1e-9, vmax=1e9)
    plt.savefig(f"{IMAGE_PATH}dm/dm_{count:03d}.png", bbox_inches="tight", dpi=300)

    pynbody.plot.image(s.g, qty="rho", width=8000, cmap="inferno", vmin=1e-8, vmax=1e7)
    plt.savefig(f"{IMAGE_PATH}gas/gas_{count:03d}.png", bbox_inches="tight", dpi=300)

    count += 1
