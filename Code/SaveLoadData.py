import pylab as py
import pickle
import numpy as np


#### load data from the file given
def loadData(namefile):
    f = open(namefile, 'rb')
    loaded_objects = []
    for i in range(11):
        loaded_objects.append(pickle.load(f))
    f.close()
    return loaded_objects[0],loaded_objects[1],loaded_objects[2],loaded_objects[3],loaded_objects[4],loaded_objects[5],loaded_objects[6],loaded_objects[7],loaded_objects[8],loaded_objects[9],loaded_objects[10]
    

#### Save data into the file given
def saveData(t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym,namefile):
    f = open(namefile, 'wb')
    for obj in [t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym]:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    f.close()


#### Fonction identique mais pour sauvegarder la representation temps/frequence. Je conseille vivement de ne pas les utiliser sauf si vous avez une raison bien precise... Ficher TRES lourd en sortie!
def saveTFR(TFR,T,F,namefile):
	f = open(namefile, 'wb')
	for obj in [TFR,T,F]:
		pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
	f.close()

def loadTFR(namefile):
	f = open(namefile, 'rb')
	loaded_objects = []
	for i in range(3):
		loaded_objects.append(pickle.load(f))
	f.close()
	return loaded_objects[0],loaded_objects[1],loaded_objects[2]

#### Test unitaire, peut etre utile si vous voulez voir comment charger ou sauvegarder des donnees ;)
if __name__ == '__main__':
	# saveData(t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym,'trans2.pckl')
	t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym = loadData('trans2.pckl')
	py.figure()
	py.plot(t,np.real(dif))    
	py.figure()
	py.scatter(x,y)  
	py.scatter(xm,ym)    
	py.show()
