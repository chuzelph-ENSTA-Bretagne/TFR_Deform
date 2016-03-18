import numpy as np
import math as m
from pyIbex import *
from vibes import vibes
import time
from SaveLoadData import *


if __name__ == '__main__':
	A=[1,1]
	B=[0,-1]
	C=[-1,0]

	vibes.beginDrawing()
	vibes.newFigure('TestCase1')
	vibes.setFigureProperties({'x': 0, 'y': 0,'width': 500, 'height': 500})
	vibes.axisLimits(-25, 25, -25, 25)
	t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym = loadData("rot1.pckl")
	t,dt,Field,inc,dif,x2,y2,xmp,ymp,xm,ym = loadData("rot2.pckl")
	time.sleep(15)
	for i in range(len(t)):

		vibes.drawCircle(x[i], y[i], 1, 'black[green]')
		vibes.drawCircle(x2[i], y2[i], 1, 'black[red]')
		vibes.drawCircle(xm, ym, 0.4, 'black[blue]')
		vibes.drawArrow([20, 15], [15, 15], 1, 'black[black]')
		vibes.drawArrow([20, 10], [15, 10], 1, 'black[black]')
		vibes.drawArrow([20, 5], [15, 5], 1, 'black[black]')
		vibes.drawArrow([20, 0], [15, 0], 1, 'black[black]')
		vibes.drawArrow([20, -5], [15, -5], 1, 'black[black]')
		time.sleep(0.001)
		vibes.clearFigure()
	vibes.endDrawing()
