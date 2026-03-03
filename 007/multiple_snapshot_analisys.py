import os
import matplotlib.pyplot as plt

SNAPSHOT_CODES = ["0005_0006_0/", 
                  "0005_0006_100/", 
                  "0005_0006_1000/", 
                  "0005_0006_2000/", 
                  "0005_0006_3000/"]
IMAGE_PATH = "plots/velocity_v0/"
FONT_SIZE = 14
os.makedirs(f"{IMAGE_PATH}", exist_ok=True)

SNAPSHOT_CONFIG = {
    "0005_0005_0/": {
        "v0": 0,
        "v_rel": 0,
        "v_shock": 0,
    },
    "0005_0005_1000/": {
        "v0": 1000,
        "v_rel": 0,
        "v_shock": 0,
    },
    "0005_0005_3000/": {
        "v0": 3000,
        "v_rel": 0,
        "v_shock": 0,
    },
    "0005_0006_0/": {
        "v0": 0,
        "v_rel": 560.6601,
        "v_shock": 1454.08,
    },
    "0005_0006_100/": {
        "v0": 100,
        "v_rel": 508.9135,
        "v_shock": 1535.23,
    },
    "0005_0006_1000/": {
        "v0": 1000,
        "v_rel": 930.8206,
        "v_shock": 1883.77,
    },
    "0005_0006_2000/": {
        "v0": 2000,
        "v_rel": 2084.6358,
        "v_shock": 2129.64,
    },
    "0005_0006_3000/": {
        "v0": 3000,
        "v_rel": 3384.9100,
        "v_shock": 2338.71,
    },
    "0005_0007_0/": {
        "v0": 0,
        "v_rel": 0,
        "v_shock": 0,
    },
    "0005_0007_1000/": {
        "v0": 1000,
        "v_rel": 0,
        "v_shock": 0,
    },
    "0005_0007_3000/": {
        "v0": 3000,
        "v_rel": 0,
        "v_shock": 0,
    },
}

plot_data = {"v0": [], "v_rel": [], "v_shock": []}
for SNAPSHOT_CODE in SNAPSHOT_CODES:
    try:
        cfg = SNAPSHOT_CONFIG[SNAPSHOT_CODE]
        plot_data["v0"].append(cfg["v0"])
        plot_data["v_rel"].append(cfg["v_rel"])
        plot_data["v_shock"].append(cfg["v_shock"])
    except KeyError:
        raise ValueError(f"Configuração não encontrada para SNAPSHOT_CODE = {SNAPSHOT_CODE}")

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(plot_data["v0"], plot_data["v_rel"], 'o--', ms = 8, label='Velocidade Relativa ($v_{rel}$)')
ax.plot(plot_data["v0"], plot_data["v_shock"], '^--', ms = 8, label='Velocidade da Frente de Choque ($v_{shock}$)')
ax.set_ylabel(r'$v_{\mathrm{rel}},\ v_{\mathrm{shock}}\ \mathrm{(km\,s^{-1})}$',
              fontsize=FONT_SIZE)
ax.set_xlabel('$v_{\mathrm{0}}\ \mathrm{(km\,s^{-1})}$', fontsize=FONT_SIZE)
ax.tick_params(axis='both', labelsize=FONT_SIZE)
ax.set_aspect('auto')
ax.legend(loc="upper left", fontsize=FONT_SIZE)
plt.tight_layout()
# plt.show()
plt.savefig(f"{IMAGE_PATH}0005_0006.png")
