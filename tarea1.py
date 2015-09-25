# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 23:44:06 2015

@author: splatt
"""

import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import integrate
plt.figure(1)
plt.clf()
longitud_de_onda=np.array([])
flujo=np.array([])
datos=np.genfromtxt('sun_AM0.dat')
for elemento in datos:
    longitud_de_onda=np.append(longitud_de_onda,elemento[0])
    flujo=np.append(flujo,elemento[1])
longitud_de_onda*=10 #Armstrong
flujo*=10000000000 #g*s-3*cm-1
plt.plot(longitud_de_onda,flujo,label='longitud de onda')
plt.xlim([longitud_de_onda[0],longitud_de_onda[len(longitud_de_onda)-1]])
plt.ylim([0,100000])
plt.xlabel('Longitud de onda [Armstrong]')
plt.ylabel('Flujo[g*s-3*cm-1]')
plt.title('Grafico Longitud de onda v/s Flujo')
plt.show()
longitud_de_onda/=10 #nm
flujo/=10000000000 #W/m²/nm

def integral_metodo_trapecio(arrx,arry,n):
    k=0
    I=0.0
    l=len(arry)
    while k<n:
        A=((arry[k*(l-1)/n])+(arry[(k+1)*(l-1)/n]))*((arrx[(k+1)*(l-1)/n])-(arrx[k*(l-1)/n]))/2.0
        I+=A
        k+=1
    return I

P2=integral_metodo_trapecio(longitud_de_onda,flujo,1696)
print "Integral P2= " + str(P2)

Luminosidad=P2*(4)*(np.pi)*((1.496*(10**11))**2)
print "Luminosidad= "+str(Luminosidad)

def funcion_de_planck(y):
    H3=6.62606957e-34**3
    K4=1.3806488e-23**4
    C2=299792458*299792458
    T4=5778**4
    a=2*math.pi*T4*K4/H3/C2
    tg3=math.tan(y)**3
    sec2=(1/math.cos(y))**2
    return a*tg3*sec2/((math.e)**(math.tan(y))-1)

H3=6.62606957e-34**3
K4=1.3806488e-23**4
C2=299792458*299792458
T4=5778**4
a=2*math.pi*T4*K4/H3/C2   
    
def integral_metodo_simpson(f,a,b,n):
    deltax=(b-a)/n
    k=0
    I=0
    while k<n:
        xk=k*deltax+a
        if k%2==0:
            I+=2*f(xk)
        else:
            I+=4*f(xk)
        k+=1
    I+=f(a)+f(b)
    I*=deltax/3
    return I

def integral_metodo_punto_medio(f,a,b,n):
    deltax=(b-a)/n
    I=0
    k=0
    while k<n:
        I+=(((a+(k+1)*deltax)-a+k*deltax))*f(((a+k*deltax)+(a+(k+1)*deltax))/2)
        k+=1
    return I
    
p="s"
while (p=="s"): 
    p3=input("n?")
    I=integral_metodo_simpson(funcion_de_planck,(math.pi/2)/p3,(math.pi/2)-(math.pi/2)/p3,p3-2)   
    I+=integral_metodo_punto_medio(funcion_de_planck,0,(math.pi/2)/p3,1)
    I+=integral_metodo_punto_medio(funcion_de_planck,(math.pi/2)-(math.pi/2)/p3,(math.pi/2),1)
    error=100-(I*100*15/((math.pi)**4)/a)
    b=str(error)
    p=raw_input("Error es "+b+"%. Cambiar n? s/n")
print "Integral P3= "+ str(I)

print "Comparación resultados de P2 y P3 (P3/P2) = "+ str(I/P2)

rs=math.sqrt(I*((1.496*(10**11))**2)/T4/(5.670373*(10**(-8))))
print "Radio solar = " + str(rs)

met1=np.trapz(flujo,longitud_de_onda)
met2=integrate.quad(funcion_de_planck,0.01,math.pi/2.01)

print "Resultado np.trapz= "+ str(met1)
print "Resultado integrate.quad= "+str(met2)
err=100-(I*100/met2[0])
print "error entre metodo scipy y mi metodo ="+str(err)