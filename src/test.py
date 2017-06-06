import numpy as np
from sklearn.datasets import fetch_mldata
from PIL import Image
import mixture_bernoulli_wrap
import sys

mnist = fetch_mldata('MNIST original', data_home=".")
zero_index = mnist["data"] <= 128
mnist["data"][zero_index] = 0
mnist["data"][np.logical_not(zero_index)] = 1
target_index = np.logical_or(np.logical_or(mnist["target"] == 2, mnist["target"] == 3), mnist["target"] == 4)
data = mnist["data"][target_index]
np.random.shuffle(data)
data = data[:600]

means = mixture_bernoulli_wrap.mixtureBernoullWrap(data, 3)
for i, mean in enumerate(means):
    mean = np.array(mean)
    pil_img = Image.fromarray(np.uint8(mean.reshape(28, 28) * 255))
    pil_img.save('result' + str(i) + '.jpg')
