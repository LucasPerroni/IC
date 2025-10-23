from PIL import Image
import matplotlib.pyplot as plt

SNAPSHOT_CODE_1 = "0005_0006_100/"
SNAPSHOT_CODE_2 = "0005_0006_1000/"
SNAPSHOT_CODE_3 = "0005_0006_2000/"
IMAGE_PATH_1 = "plots/" + SNAPSHOT_CODE_1
IMAGE_PATH_2 = "plots/" + SNAPSHOT_CODE_2
IMAGE_PATH_3 = "plots/" + SNAPSHOT_CODE_3

# DIST X TIME
pos_time = [f"{IMAGE_PATH_1}dist.png", f"{IMAGE_PATH_2}dist.png", f"{IMAGE_PATH_3}dist.png"]
imgs = [Image.open(arq) for arq in pos_time]

largura_total = sum(img.width for img in imgs)
altura = max(img.height for img in imgs)
nova = Image.new("RGB", (largura_total, altura), (255, 255, 255))

x_offset = 0
for im in imgs:
    nova.paste(im, (x_offset, 0))
    x_offset += im.width
nova.save("plots/final/dist.png")

# VELOCITY X POS/TIME
t_r = [
    f"{IMAGE_PATH_1}vel-pos.png", f"{IMAGE_PATH_2}vel-pos.png", f"{IMAGE_PATH_3}vel-pos.png",
    f"{IMAGE_PATH_1}vel-time.png", f"{IMAGE_PATH_2}vel-time.png", f"{IMAGE_PATH_3}vel-time.png"
]

imgs = [Image.open(arq) for arq in t_r]
w, h = imgs[0].size  

cols, rows = 3, 2
nova = Image.new("RGB", (cols * w, rows * h), (255, 255, 255))

for i, im in enumerate(imgs):
    col = i % cols
    row = i // cols
    nova.paste(im, (col * w, row * h))
nova.save("plots/final/vel-pos_time.png")
