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
    r=mth.sqrt(tau**2 + ((sx-sy)/2)**2)
    s_p=(sx+sy)/2
    c=(s_p,0)
    m = -2*tau/sx
    X = np.linspace(-r+s_p,r+s_p)
    Y = m*X - m*s_p
    plt.figure(dpi = 500)
    figure, axes = plt.subplots()
    draw_circle = plt.Circle(c, r,fill=False)
    plt.plot(X,Y)
    plt.xlim(-2.5*r,2.5*r)
    plt.ylim(-2.5*r,2.5*r)
    plt.hlines(0,-2.5*r,2.5*r,'k')
    plt.vlines(0,-2.5*r,2.5*r,'k')
    axes.set_aspect(1)
    axes.add_artist(draw_circle)
    plt.title('Circulo Mohr')
    plt.grid(axis='both')
    plt.show()
    return None
    


