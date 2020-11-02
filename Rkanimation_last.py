# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 17:55:07 2020

@author: Lenovo
"""


import numpy as np
from copy import copy
from func import function_evaluation
from func import energy_calculation


"""
This is a code for general explicit RK method. Classical RK method 
and Explicit Euler method are sub-categories of RK method.
We have to change alpha,gamma,beta arrays accordingly. together alpha,beta 
and gamma are called the butcher table
"""

#defining parameters of RK method
alpha=np.array([0,1/2,1/2,1]);
gamma=np.array([1/6,1/3,1/3,1/6]);
beta=np.array([[0,0,0,0],[0.5,0,0,0],[0,0.5,0,0],[0,0,1,0]])

"""
t0=starting time, T=maximum time, h=step size, y0=inital condition array
tvec stores value of time, yvec stores yj at each time step. yj is a 4x1 array
with the first entry theta1,second entry theta2, third entry p1
 and the fourth entry q2
As this is an explicit method, we don't need any Newton Raphson iteration scheme 
"""
#setting inital parameters
t0=0;
T=10;
h=0.01;
y0=np.array([1.04,1.04,0,0]).reshape(-1,1);

#setting parameters of pendulum
m1=2;m2=1;l1=1;l2=1;

num_steps=round((T-t0)/h);

#values of the four variables will be stored in yvec for each time step and the corresponding time step will be stored in tvec

yvec = np.empty((4,0))
tvec=np.empty((1,0))


yj=y0.copy()
tj=copy(t0)

"""
getting the values of RK method from the butcher table
m indicates the order of the RK method
num_eqts indicate the number of equations in the RHS of the ODE.
in our case we have four equations in the RHSof the ODE. 

"""
m=np.size(alpha)
num_eqts=np.size(y0)


print("Calculating...")
counter=1;
while counter<=num_steps:
    yvec=np.hstack((yvec,yj))
    tvec=np.append(tvec,tj)
    
    q1=yj[0,0]
    q2=yj[1,0]
    p1=yj[2,0]
    p2=yj[3,0]
       
    x1,x2,x3,x4=function_evaluation(q1,q2,p1,p2,m1,m2,l1,l2)
  
    """
    k matrix is used in RK method calculation.
    """
    
    #initializing k matrix
    k=np.zeros(shape=(num_eqts,m))
    #startting the loop to calculate beta*k term
    i=0
    while i<m:
        beta_k=np.zeros((num_eqts,1))
        l=0
        sum_term=0
        while l<i:
            sum_term=beta[i,l]*k[:,l]
            beta_k=beta_k+np.array(sum_term).reshape(-1,1)
            l+=1
        #'l' loop ends here   
        
        
    #updating the values of k matrix using the given formula
        h_beta_k=yj+h*beta_k
        q1=h_beta_k[0,0]
        q2=h_beta_k[1,0]
        p1=h_beta_k[2,0]
        p2=h_beta_k[3,0]
        
        x1,x2,x3,x4=function_evaluation(q1,q2,p1,p2,m1,m2,l1,l2)
        
        k[0,i]=x1
        k[1,i]=x2
        k[2,i]=x3
        k[3,i]=x4
        
        
        i+=1
    #"i" loop ends here
    
    
    gamma_k=np.beta_k=np.zeros(shape=(num_eqts,1))
    l=0;
    sum_term=0
    
    while l<m:
        sum_term=gamma[l]*k[:,l]
        gamma_k=gamma_k+np.array(sum_term).reshape(-1,1)
        l+=1
    #last inner loop ends here
    
    yj=yj+h*gamma_k;
    tj=tj+h;
    counter+=1

"""
Calculating Energy, Velocity and Acceleration
"""

pot_ener_vec=np.zeros(num_steps)
kin_ener_vec=np.empty(num_steps)
vel_vec = np.zeros((2,num_steps))
acc_vec=np.zeros((2,num_steps))

i=0
while i<num_steps:
    
    yj=yvec[:,i]
    q1=yj[0]
    q2=yj[1]
    p1=yj[2]
    p2=yj[3]
   
    x1,x2,x3,x4=function_evaluation(q1,q2,p1,p2,m1,m2,l1,l2)
    P,T=energy_calculation(q1,q2,p1,p2,m1,m2,l1,l2,x1,x2,x3,x4)
    pot_ener_vec[i]=P
    kin_ener_vec[i]=T
    vel_vec[0,i]=x1
    vel_vec[1,i]=x2
    acc_vec[0,i]=x3/m1
    acc_vec[1,i]=x4/m2
    
    i+=1

total_energy=pot_ener_vec+kin_ener_vec  


# """
# PLOTTING
# """ 
# print("plotting")
# import matplotlib.pyplot as plt
# dev_x=tvec
# plt.figure(0)
# plt.plot(dev_x,total_energy)
# plt.plot(dev_x,kin_ener_vec)
# plt.plot(dev_x,pot_ener_vec)
# plt.xlabel("Time")
# plt.ylabel("Energy")
# plt.legend(["Total Energy","Kinetic Energy","Potential Energy"])
# plt.title("Conservation of Energy Over Time Using Classical RK")

# plt.figure(1)
# plt.plot(dev_x,yvec[0,:])
# plt.plot(dev_x,yvec[1,:])
# plt.xlabel("Time")
# plt.ylabel("Position")
# plt.legend(["Mass 1","Mass 2"])
# plt.title("Position of Mass1 and Mass 2 Using Classical RK")

# plt.figure(2)
# plt.plot(dev_x,vel_vec[0,:])
# plt.plot(dev_x,vel_vec[1,:])
# plt.xlabel("Time")
# plt.legend(["Mass 1","Mass 2"])
# plt.title("Velocity of Mass1 and Mass 2 Using Classical RK")

# plt.figure(3)
# plt.plot(dev_x,acc_vec[0,:])
# plt.plot(dev_x,acc_vec[1,:])
# plt.xlabel("Time")
# plt.legend(["Mass 1","Mass 2"])
# plt.title("Acceleration of Mass1 and Mass 2 Using Classical RK")

"""

"""
import pygame
import sys
from pygame.locals import *
from math import sin, cos, pi
from numpy.linalg import inv
from funcanimation import y_update


def update(a1, a2):
      scale = 100
      offset = (400, 50)
      global m1;m2;l1;l2
      x1 = l1*scale*sin(a1)+offset[0]
      y1 = l1*scale*cos(a1)+offset[1]
      x2 = x1+l2*scale*sin(a2)
      y2 = y1+l2*scale*cos(a2)

      return (x1, y1),(x2, y2)
 
    
def render(point1, point2):
      scale = 10
      offset = (400, 50)
      x1, y1,  = int(point1[0]), int(point1[1])
      x2, y2,  = int(point2[0]), int(point2[1])

      if prev_point:
        xp, yp = prev_point[0], prev_point[1]
        pygame.draw.line(trace, LT_BLUE, (xp, yp), (x2, y2), 3)

      screen.fill(WHITE)    
      screen.blit(trace, (0,0))

      pygame.draw.line(screen, BLACK, offset, (x1,y1), 5)
      pygame.draw.line(screen, BLACK, (x1,y1), (x2,y2), 5)
      pygame.draw.circle(screen, BLACK, offset, 8)
      pygame.draw.circle(screen, MAGENTA, (x1, y1), int(m1*scale))
      pygame.draw.circle(screen, BLUEGREEN, (x2, y2), int(m2*scale))

      return (x2, y2)

wide, height = 800, 300
WHITE = (255,255,255)
MAGENTA =(255,0,230)
MOON = (235,245,100)
BLUEGREEN =(0,255,170)
BLACK = (0,0,0)
RED = (255,0,0)
BLUE = (0,0,255)
LT_BLUE = (230,230,255)
offset = (400, 50)

screen = pygame.display.set_mode((wide,height))
screen.fill(WHITE)
trace = screen.copy()
pygame.display.update()
clock = pygame.time.Clock()

prev_point = None
t = 0
print(h)
dt= h
y = y0

pygame.font.init()
myfont = pygame.font.SysFont('Calibri', 25)



while t<=10:
    
  
  for event in pygame.event.get():
  		if event.type == pygame.QUIT:
   			sys.exit()
  
 
  point1, point2 = update(y[0], y[1])
  prev_point = render(point1, point2)
 
  
 
  t += dt
  #print(t)
  index = np.where(tvec==t)
  print(index)
  if len(index[0])!=1:
      print("done")
      break
      
  else:
      i= int(index[0])
      print(type(index[0]))
 
      y = yvec[:,i]
      
      
      time_string = 'Time: {} seconds'.format(round(t,1))
      text = myfont.render(time_string, False, (0, 0, 0))
      screen.blit(text, (10,10))
      
      
      clock.tick(60)
      pygame.display.update()
  
pygame.display.quit()
pygame.quit()
sys.exit()         