from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt

SNAPSHOT_CODE_1 = "0005_0006_100/"
SNAPSHOT_CODE_2 = "0005_0006_1000/"
SNAPSHOT_CODE_3 = "0005_0006_2000/"
IMAGE_PATH_1 = "plot/" + SNAPSHOT_CODE_1
IMAGE_PATH_2 = "plot/" + SNAPSHOT_CODE_2
IMAGE_PATH_3 = "plot/" + SNAPSHOT_CODE_3

# POSITION X TIME
pos_time = [f"{IMAGE_PATH_1}pos-time.png", f"{IMAGE_PATH_2}pos-time.png", f"{IMAGE_PATH_3}pos-time.png"]
imgs = [Image.open(arq) for arq in pos_time]

largura_total = sum(img.width for img in imgs)
altura = max(img.height for img in imgs)
nova = Image.new("RGB", (largura_total, altura), (255, 255, 255))

x_offset = 0
for im in imgs:
    nova.paste(im, (x_offset, 0))
    x_offset += im.width
nova.save("plot/final/pos_time.png")

# TEMPERATURE X RADIUS
t_r = [
    f"{IMAGE_PATH_1}t-r_090.png", f"{IMAGE_PATH_2}t-r_060.png", f"{IMAGE_PATH_3}t-r_040.png",
    f"{IMAGE_PATH_1}t-r_094.png", f"{IMAGE_PATH_2}t-r_064.png", f"{IMAGE_PATH_3}t-r_044.png",
    f"{IMAGE_PATH_1}t-r_098.png", f"{IMAGE_PATH_2}t-r_068.png", f"{IMAGE_PATH_3}t-r_048.png"
]

imgs = [Image.open(arq) for arq in t_r]
w, h = imgs[0].size  

cols, rows = 3, 3
nova = Image.new("RGB", (cols * w, rows * h), (255, 255, 255))

for i, im in enumerate(imgs):
    col = i % cols
    row = i // cols
    nova.paste(im, (col * w, row * h))
nova.save("plot/final/t-r.png")
