import numpy as np
import math as m
from pyIbex import *
from vibes import vibes
import time


if __name__ == '__main__':
	A=[1,1]
	B=[0,-1]
	C=[-1,0]

	vibes.beginDrawing()
	vibes.setFigureProperties({'x': 0, 'y': 0,'width': 500, 'height': 500})
	vibes.axisLimits(-25, 25, -25, 25)

	for t in range(500):

		vibes.drawCircle(np.cos(2*np.pi*(t/360)), np.sin(2*np.pi*(t/360)), 0.4, 'black[green]')
		vibes.drawCircle(-np.sin(2*np.pi*(t/360)), np.cos(2*np.pi*(t/360)), 0.4, 'black[green]')
		vibes.drawCircle(-1+t*0.01, -1+t*0.01, 0.4, 'black[green]')
		time.sleep(0.1)
		vibes.clearFigure()
	vibes.endDrawing()
