#!/usr/bin/env python3
import sys
sys.path.append('/usr/local/lib/python3.8.10/site-packages')
import numpy as np
from array import *
from matplotlib.colors import colorConverter
import matplotlib as mpl
import matplotlib.pyplot as plt  
from PIL import Image

im= Image.open('Images/renfe_lineas.png','r') #or metro_lineas.png

im = im.resize((500, 500))
width=im.size[0] 
height=im.size[1] 
tl = 250
num = 0
for i in range(0, width): 
    for j in range(0, height):
        data=im.getpixel((i,j))
        if (int(data[0])>tl and int(data[1])>tl and int(data[2])>tl):
            num=num+1
pp=num/(500*500)
print(pp)
