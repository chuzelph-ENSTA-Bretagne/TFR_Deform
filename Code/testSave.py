
import pylab as py
import pickle
import numpy as np


def loadData(namefile):
    f = open(namefile, 'rb')
    loaded_objects = []
    for i in range(4):
        loaded_objects.append(pickle.load(f))
    f.close()
    return loaded_objects[0],loaded_objects[1],loaded_objects[2],loaded_objects[3]
    
def saveData(t,listData,x,y,namefile):
    f = open(namefile, 'wb')
    for obj in [t , listData,x,y]:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    f.close()


if __name__ == '__main__':
    t,List2,x,y = loadData('Tran2.pckl')
    py.figure()
    py.plot(t,np.real(List2))    
    py.figure()
    py.scatter(x,y)    
    py.show()