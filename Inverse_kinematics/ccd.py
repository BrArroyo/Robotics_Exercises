#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Robótica Computacional - 
# Grado en Ingeniería Informática (Cuarto)
# Práctica: Resolución de la cinemática inversa mediante CCD
#           (Cyclic Coordinate Descent).

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import colorsys as cs

# ******************************************************************************
# Declaración de funciones

def muestra_origenes(O,final=0):
  # Muestra los orígenes de coordenadas para cada articulación
  print('Origenes de coordenadas:')
  for i in range(len(O)):
    print('(O'+str(i)+')0\t= '+str([round(j,3) for j in O[i]]))
  if final:
    print('E.Final = '+str([round(j,3) for j in final]))

def muestra_robot(O,obj):
  # Muestra el robot graficamente
  plt.figure(1)
  plt.xlim(-L,L)
  plt.ylim(-L,L)
  T = [np.array(o).T.tolist() for o in O]
  for i in range(len(T)):
    plt.plot(T[i][0], T[i][1], '-o', color=cs.hsv_to_rgb(i/float(len(T)),1,1))
  plt.plot(obj[0], obj[1], '*')
  plt.show()
  raw_input()
  plt.clf()

def matriz_T(d,th,a,al):
  # Calcula la matriz T (ángulos de entrada en grados)
  
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]

def cin_dir(th,a):
  #Sea 'th' el vector de thetas
  #Sea 'a'  el vector de longitudes
  T = np.identity(4)
  o = [[0,0]]
  for i in range(len(th)):
    T = np.dot(T,matriz_T(0,th[i],a[i],0))
    tmp=np.dot(T,[0,0,0,1])
    o.append([tmp[0],tmp[1]])
  return o

# ******************************************************************************
# Cálculo de la cinemática inversa de forma iterativa por el método CCD

# valores articulares arbitrarios para la cinemática directa inicial
th=[0.,0.,0.,0.]
a =[5.,5.,5.,2.]
prismatica = [False, False, False, True]
L = sum(a) # variable para representación gráfica
EPSILON = .01

# Limites de los articulaciones
lim_max = [pi/2, pi/2, pi/2, 10]
lim_min = [-pi/2, -pi/2, -pi/2, 0]

plt.ion() # modo interactivo

# introducción del punto para la cinemática inversa
if len(sys.argv) != 3:
  sys.exit("python " + sys.argv[0] + " x y")
objetivo=[float(i) for i in sys.argv[1:]]

O=range(len(th)+1) # Reservamos estructura en memoria
O[0]=cin_dir(th,a) # Calculamos la posicion inicial
print ("- Posicion inicial:")
muestra_origenes(O[0])

dist = float("inf")
prev = 0.
iteracion = 1
while (dist > EPSILON and abs(prev-dist) > EPSILON/100.):
  prev = dist
  # Para cada combinación de articulaciones:
  for i in range(len(th)):
    pos = len(th) - i -1
    # cálculo de la cinemática inversa:
    #if es una articulación de revolución 
    if not prismatica[pos]:
      # Calculamos el angulo 
      v1 = [O[i][len(th)][0]-O[i][pos][0],O[i][len(th)][1] - O[i][pos][1]]
      v2 = [objetivo[0] - O[i][pos][0],objetivo[1] - O[i][pos][1]]
      alfa1 = atan2(v1[1],v1[0])
      alfa2 = atan2(v2[1],v2[0])
      
      # Actualizamos el valor de la articulación
      th[pos] = th[pos] + alfa2 - alfa1  
    
      # corregimos los angulos con los límites
      if th[pos] > lim_max[pos]:
        th[pos] = lim_max[pos]
      if th[pos] < lim_min[pos]:
        th[pos] = lim_min[pos]  
      
    #else es una articulación prismática
    else:
      # Calculamos la suma de las articulaciones
      omega = sum(th[0:pos])
      
      # Calculamos los vectores correspondientes y calculamos el producto escalar
      vector_despesplazamiento = [cos(omega),sin(omega)]
      vector_acercamiento = np.subtract(objetivo, [O[i][pos][0], O[i][pos][1]])
      distancia = np.dot(vector_despesplazamiento, vector_acercamiento)
      
      # Actualizamos el valor de la articulación
      a[pos] += distancia
      
      # corregimos las distancia con los limites
      if a[pos] > lim_max[pos]:
        a[pos] = lim_max[pos]
      if a[pos] < lim_min[pos]:
        a[pos] = lim_min[pos] 
      
    O[i+1] = cin_dir(th,a)

  dist = np.linalg.norm(np.subtract(objetivo,O[-1][-1]))
  print ("\n- Iteracion " + str(iteracion) + ':')
  muestra_origenes(O[-1])
  muestra_robot(O,objetivo)
  print ("Distancia al objetivo = " + str(round(dist,5)))
  iteracion+=1
  O[0]=O[-1]

if dist <= EPSILON:
  print ("\n" + str(iteracion) + " iteraciones para converger.")
else:
  print ("\nNo hay convergencia tras " + str(iteracion) + " iteraciones.")
print ("- Umbral de convergencia epsilon: " + str(EPSILON))
print ("- Distancia al objetivo:          " + str(round(dist,5)))
print ("- Valores finales de las articulaciones:")
for i in range(len(th)):
  print ("  theta" + str(i+1) + " = " + str(round(th[i],3)))
for i in range(len(th)):
  print ("  L" + str(i+1) + "     = " + str(round(a[i],3)))
