
import pylab as py
import pickle
import numpy as np
import tftb
import matplotlib.pyplot as plt   
import myspectrogram2


def loadData(namefile):
    f = open(namefile, 'rb')
    loaded_objects = []
    for i in range(11):
        loaded_objects.append(pickle.load(f))
    f.close()
    return loaded_objects[0],loaded_objects[1],loaded_objects[2],loaded_objects[3],loaded_objects[4],loaded_objects[5],loaded_objects[6],loaded_objects[7],loaded_objects[8],loaded_objects[9],loaded_objects[10]
    
def saveData(t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym,namefile):
    f = open(namefile, 'wb')
    for obj in [t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym]:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    f.close()


if __name__ == '__main__':
	t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym = loadData('trans2.pckl')
	py.figure()
	py.plot(t,np.real(dif))    
	py.figure()
	py.scatter(x,y)    
	py.show()
	tfr = tftb.processing.ShortTimeFourierTransform(np.real(dif),n_fbins = 4000,fwindow = np.hamming(2000))		#n_fbins
	tfr.run()
	tfr.plot(kind='cmap', show_tf=True)



	#X , T , F = myspectrogram2.myspectrogram(np.real(dif), 1000, 2, 1000, 10000)
	#plt.figure()
	#plt.imshow(abs(X), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
	#plt.xlabel('Time')
	#plt.ylabel('Frequency')
	#plt.colorbar()
	#plt.show()

	#wvd = tftb.processing.cohen.WignerVilleDistribution(np.real(dif))
	#wvd.tfr
	#wvd.run()
	#wvd.plot(kind='cmap', show_tf=True) 
	#plt.figure()
	#plt.imshow(abs(wvd.tfr), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
	#plt.xlabel('Time')
	#plt.ylabel('Frequency')
	#plt.colorbar()
	#plt.show()
