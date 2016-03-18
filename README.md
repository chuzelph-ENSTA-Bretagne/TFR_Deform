# ApplicationSystem54
Résultat de mon application système de l'UV 5.4 / Resultat of the system application UV 5.4

LawFrequency.py : Contains 3 functions

- splineFit(lawFound,t,smooth) : Realize the interpolation of lawFound, according the vector t and a smoothing coefficient smooth.

- seamCarving(tfr,t,f) : use the seam carving algorithm of the time frequency representation tfr, according the vector t and f in order to return the law found. Use lawFrequency(tfr,t,f).

- lawFrequency(tfr,t,f) : Implement the seam carving algorithm.

ObjectDiff.py: Script which realize the generation of the signal. This class create an objet which can generate and then save a signal according to specific parameters that have to be given by the oprator. 

GenerationSignal.py : This class implement many functions but the most important is Field_scat(k,r0,r,alpha,theta,N). It generate the field scattered by the cylinder according to the value of the radius r0, the distance r and the parameters alpha, theta, k and N which is the number of terms we will keep in the serie.

TFR_Signal.py : This script show many functions which can be use in order to realize the time frequency representation of a signal.

# Version 1.3

All the report is over. The code can generate and analyse a signal. In order to easily use the code, two script are created and can be lunch with the command :

sh "name of the script"

Those script will first generate some signal which are representative of the kind of signal we have to analyse and secondly will expose how the interpolation works.

# Version 1.2

This version implemente all what is needed in order to find the most significative frequency law inside a signal and find an interpolation. We can easily see the impact of the movement of an object on the electronic signature (Doppler effect).

# Version 1.1

Because of an issue with github, I had to create a new github folder. In this version, I implemented an script which can easily generate the filed generate by a cylinder due to an electromagnetic wave. 

# Version 1.0

This version is the result of my preparation during the month of october to january. I have implemented my own function for the time frequency representation. Nevertheless, I recently found a toolbox similar to the Flandrin toolbox in Matlab. I deceided to use its function for this project.
