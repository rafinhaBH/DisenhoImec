''' Código para resolver la parte 2 del laboratorio 1 de disismec Modelo de pared delgada con cilindro abierto'''
import matplotlib.pyplot as plt
import numpy as np
r = 36.5/1000 #m
t = 3e-6 #m
v = 0.3
E = 72e+9
import stressTransformations as st 

### Para P1 = 10 Bar
p1 = 10e+5
ec1 = -60e-6 #Para quedar en m/m
eb1 = 60e-6
ea1 = 173e-6
# Cálculo de deformaciones principales
e1_1 = (ea1+eb1)/2 + 0.5*np.sqrt(2)*np.sqrt( ((ea1 - eb1)**2) + ((eb1-ec1)**2) )
e2_1 = (ea1+eb1)/2 - 0.5*np.sqrt(2)*np.sqrt( ((ea1 - eb1)**2) + ((eb1-ec1)**2) )
# Cálculos de esfuerzos principales
S1_1 = ((E/(1-v*v))*(e1_1+v*e2_1))*1e-6
S1_2 = ((E/(1-v*v))*(e2_1+v*e1_1))*1e-6
## Circulos de mohr
radio1 = (S1_1 - S1_2)/2
c1 = (S1_1+S1_2)/2
circle1 = plt.Circle((c1, 0),radio1, fill=False)
fig1, ax = plt.subplots()
fig1.dpi=1000
plt.xlim(0,3*radio1)
plt.ylim(-1.5*radio1,1.5*radio1)
plt.grid(linestyle = '--')
ax.add_artist(circle1)
ax.set_aspect(1)
plt.title('Circulo de Mohr [MPa] para P1 = 10bar con deformaciones')
plt.show()
print(c1,0)

### Para p2 = 20 Bar
p2 = 20e+5
ec2 = -124e-6
eb2 =  117e-6
ea2 =  352e-6
# Cálculo de deformaciones principales
e1_2 = (ea2+eb2)/2 + 0.5*np.sqrt(2)*np.sqrt( ((ea2 - eb2)**2) + ((eb2-ec2)**2) )
e2_2 = (ea2+eb2)/2 - 0.5*np.sqrt(2)*np.sqrt( ((ea2 - eb2)**2) + ((eb2-ec2)**2) )
# Cálculos de esfuerzos principales
S2_1 = ((E/(1-v*v))*(e1_2+v*e2_2))*1e-6
S2_2 = ((E/(1-v*v))*(e2_2+v*e1_2))*1e-6
## Circulos de mohr
radio2 = (S2_1 - S2_2)/2
c2 = (S2_1+S2_2)/2
circle1 = plt.Circle((c2, 0),radio2, fill=False)
fig2, ax = plt.subplots()
fig2.dpi=1000
plt.xlim(0,3*radio2)
plt.ylim(-1.5*radio2,1.5*radio2)
plt.grid(linestyle = '--')
ax.add_artist(circle1)
ax.set_aspect(1)
plt.title('Circulo de Mohr [MPa] para P2 = 20bar con deformaciones')
plt.show()
print(c2,0)
