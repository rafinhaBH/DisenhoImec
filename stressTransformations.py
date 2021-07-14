#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This is a module of the calculations of stress transformations and Mohr's Circle
"""
Created on Thu Nov 19 17:29:26 2020

@author: rafaelbeltranhernandez

Bogotá, Colombia 

Contact: rj.beltran@uniandes.edu.co
"""

import math as mth

import matplotlib.pyplot as plt

import numpy as np

# New components of Stress

def Sigx (sx,sy,tau,th):
    """returns the \sigma_x' component of stress """
    return (  (sx+sy)/2 + mth.cos(2*mth.radians(th))*(sx-sy)/2 + tau*mth.sin(2*mth.radians(th)) )

def Sigy (sx,sy,tau,th):
    """returns the \sigma_y' component of stress """
    return (  (sx+sy)/2 - mth.cos(2*mth.radians(th))*(sx-sy)/2 - tau*mth.sin(2*mth.radians(th)) )

def TauP (sx,sy,tau,th):
    """Returns the \tau' component of srtess"""
    return ( mth.sin(2*mth.radians(th))*(-sx-sy)/2 + tau*mth.cos(2*mth.radians(th))  )

# Maxes and mins
    
def SigmaPr(sx,sy,tau):
    
    """ Magnitude of the Principal Stress """
    max=(sx+sy)/2 + mth.sqrt( tau**2 + ((sx-sy)/2)**2 )
    min=(sx+sy)/2 - mth.sqrt( tau**2 + ((sx-sy)/2)**2 )
    return (max,min)

def thPr(sx,sy,tau):
    """ Direction of the Principal Stress """
    theta=0.5*mth.atan(2*tau/(sx-sy))
    return (mth.degrees(theta))

def TauMax(sx,sy,tau):
    """ Maximum Shear Stress at 45° """
    return (mth.sqrt( tau**2 + ((sx-sy)/2)**2 ))

def mohr(sx,sy,tau):
    """ Plots the mohr's circle given the constants """
    import matplotlib.figure as fig
    r=mth.sqrt(tau**2 + ((sx-sy)/2)**2)
    s_p=(sx+sy)/2
    c=(s_p,0)
    m = -2*tau/sx
    X = np.linspace(-r+s_p,r+s_p)
    Y = m*X - m*s_p
    figure, axes = plt.subplots()
    figure.dpi=500
    draw_circle = plt.Circle(c, r,fill=False)
    plt.plot(X,Y)
    plt.xlim(-3*r,3*r)
    plt.ylim(-2.5*r,2.5*r)
    plt.hlines(0,-2.5*r,2.5*r,'k')
    plt.vlines(0,-2.5*r,2.5*r,'k')
    axes.set_aspect(1)
    axes.add_artist(draw_circle)
    plt.title('Circulo Mohr')
    plt.grid(axis='both')
    plt.show()
    return None
    

# Rotating Rings

def ts_rotating(rho,nu,w,ri,ro,r):
    ''' Returns the tangential stress of a rotating ring '''
    f1 = rho * ((3+nu)/8) * (w**2)
    f2 = (ri**2) + (ro**2) + (ri*ri*ro*ro/(r**2)) - (r**2)*(1+3*nu)/(3+nu)
    return (f1*f2)

def rs_rotating(rho,nu,w,ri,ro,r):
    ''' Returns the radial stress of a rotating ring'''
    f1 = rho * ((3+nu)/8) * (w**2)
    f2 = (ri**2) + (ro**2) - (ri*ri*ro*ro/(r**2)) - (r**2)
    return(f1*f2)

# Press amd shrink fits

def p_fit_2mats(E,nu,d,ri,ro,r):
    ''' Returns the pressure value of a fit of two different materials. Enter E and nu as tuple'''
    e0 = (1/E[0])*(nu[0]+ ( (ro**2) + (r**2)  )/((ro**2) - (r**2)) )
    e1 = (1/E[1])*(nu[1]+ ( (ri**2) + (r**2)  )/((r**2) - (ri**2)) )
    frac = ( d /(e0+e1)  )
    return (frac/r)

def p_fit(E,nu,d,ri,ro,r):
    ''' Returns pressure value of a fit for the same material'''
    f1 = E*d/(2*(r**3))
    f2 = (  (ro*ro-r*r)*(r*r-ri*ri)/(ro*ro-ri*ri)    )
    return (f1*f2)

# Contact stresses 

def p_max_sphere(E,nu,d1,d2,f):
    ''' Returns the maximum value of between two spheres'''
    frac1 = 3*f/8
    n0 = (1/E[0])*(1-(nu[0]**2))
    n1 = (1/E[1])*(1-(nu[1]**2))
    den = 1/d1 + 1/d2
    frac2 = (n0+n1)/den
    adentro = frac1 * frac2
    a = adentro**(1/3)
    pmax = 3*f/(2*np.pi* (a**2))
    return (pmax)


def stresses_contact_spheres(E,nu,d1,d2,f):
    ''' Returns the functions of stresses when two spheres have contact'''
    frac1 = 3*f/8
    n0 = (1/E[0])*(1-(nu[0]**2))
    n1 = (1/E[1])*(1-(nu[1]**2))
    den = 1/d1 + 1/d2
    frac2 = (n0+n1)/den
    adentro = frac1 * frac2
    a = adentro**(1/3)
    c = -p_max_sphere(E,nu,d1,d2,f)
    sx0 = lambda z: c* ( ( (1-abs(z/a)*np.arctan( a/z )) * (1+nu[0]) )   -    (  1/ (2* (1+z*z/(a*a)) )   )   )
    sx1 = lambda z: c* ( ( (1-abs(z/a)*np.arctan( a/z )) * (1+nu[1]) )   -    (  1/ (2* (1+z*z/(a*a)) )   )   )
    sz = lambda z: ( c / (1+z*z/(a*a)) )
    return (sx0,sx1,sz)

def p_max_cylinder(E,nu,d1,d2,f,l):
    frac1 = 2*f/(np.pi*l)
    n0 = (1/E[0])*(1-(nu[0]**2))
    n1 = (1/E[1])*(1-(nu[1]**2))
    den = 1/d1 + 1/d2
    frac2 = (n0+n1)/den
    adentro = frac1 * frac2
    b = adentro**(1/2)
    pmax = 2*f/(np.pi*b*l)
    return (pmax)

def stresses_contact_cylinders(E,nu,d1,d2,f,l):
    frac1 = 2*f/(np.pi*l)
    n0 = (1/E[0])*(1-(nu[0]**2))
    n1 = (1/E[1])*(1-(nu[1]**2))
    den = 1/d1 + 1/d2
    frac2 = (n0+n1)/den
    adentro = frac1 * frac2
    b = adentro**(1/2)
    c = -p_max_cylinder(E,nu,d1,d2,f,l)

    sx = lambda z: c*2*nu[0]* (np.sqrt (1+z*z/(b**2)) -abs(z/b) ) 
    sy = lambda z: c* (( (1+2*z*z/(b**2))/(np.sqrt (1+z*z/(b**2))) ) - 2*abs(z/b) )
    sz = lambda z: (c/np.sqrt(1+z*z/(b**2)))
    return(sx,sy,sz)

def vonMisses(sx,sy,tau):
    '''Returns the von misses stress for the plane state of stress'''
    ans = np.sqrt((sx**2) - (sx*sy) + sy**2 + 3*(tau**2))
    return (ans)


