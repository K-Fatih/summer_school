from math import sin,cos
import numpy as np


"""
function_evaluation evaluates the values of the given four hamiltonian
 equations 
"""


def function_evaluation(q1,q2,p1,p2,m1,m2,l1,l2):
    
    global g
    g=9.81
  
    x1=(-l1*p2*cos(q1 - q2) + l2*p1)/(l1**2*l2*(m1 + m2*sin(q1 - q2)**2))   
    x2=(l1*p2*(m1 + m2) - l2*m2*p1*cos(q1 - q2))/(l1*l2**2*m2*(m1 + m2*sin(q1 - q2)**2))
    
    A1 = (p1*p2*sin(q1-q2))/(l1*l2*(m1+m2*(sin(q1-q2))**2))
    
    A2 = (p1**2*m2*l2**2-2*p1*p2*m2*l1*l2*cos(q1-q2)+ p2**2*(m1+m2)*l1**2)*(sin(2*(q1-q2)))/(2*l1**2*l2**2*(m1+m2*(sin(q1-q2))**2)**2)
    
    x3 = -(m1+m2)*g*l1*sin(q1) - A1 + A2
    
    x4 = -m2*g*l2*sin(q2) + A1 - A2  
    
    return x1,x2,x3,x4

"""
der_func evaluates the inverse of tangent of the newton-raphson scheme
I calculated derivatives of the hamiltonian equations using sympy library and 
copy-pasted the results here.

d11 is the derivative of theta_1_dot wrt theta 1
d12 is the derivative of theta_1_dot wrt theta 2
and so on
"""
    
def der_func(q1,q2,p1,p2,m1,m2,l1,l2):
    
    #d1,d2,d3,d4 are the derivatives of our equations
    g=9.81
    h=0.0005
    
    d11=p2*sin(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) - 2*m2*(-l1*p2*cos(q1 - q2) + l2*p1)*sin(q1 - q2)*cos(q1 - q2)/(l1**2*l2*(m1 + m2*sin(q1 - q2)**2)**2)
    d12=-p2*sin(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) + 2*m2*(-l1*p2*cos(q1 - q2) + l2*p1)*sin(q1 - q2)*cos(q1 - q2)/(l1**2*l2*(m1 + m2*sin(q1 - q2)**2)**2)
    d13=1/(l1**2*(m1 + m2*sin(q1 - q2)**2))
    d14=-cos(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2))
    
    d21=p1*sin(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) - 2*(l1*p2*(m1 + m2) - l2*m2*p1*cos(q1 - q2))*sin(q1 - q2)*cos(q1 - q2)/(l1*l2**2*(m1 + m2*sin(q1 - q2)**2)**2)
    d22=-p1*sin(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) + 2*(l1*p2*(m1 + m2) - l2*m2*p1*cos(q1 - q2))*sin(q1 - q2)*cos(q1 - q2)/(l1*l2**2*(m1 + m2*sin(q1 - q2)**2)**2)
    d23=-cos(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2))
    d24=(m1 + m2)/(l2**2*m2*(m1 + m2*sin(q1 - q2)**2))
    
    d31=g*l1*(-m1 - m2)*cos(q1) + 2*m2*p1*p2*sin(q1 - q2)**2*cos(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)**2) + m2*p1*p2*sin(q1 - q2)*sin(2*q1 - 2*q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)**2) - p1*p2*cos(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) - 2*m2*(l1**2*p2**2*(m1 + m2) - 2*l1*l2*m2*p1*p2*cos(q1 - q2) + l2**2*m2*p1**2)*sin(q1 - q2)*sin(2*q1 - 2*q2)*cos(q1 - q2)/(l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**3) + (l1**2*p2**2*(m1 + m2) - 2*l1*l2*m2*p1*p2*cos(q1 - q2) + l2**2*m2*p1**2)*cos(2*q1 - 2*q2)/(l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**2)
    d32=-2*m2*p1*p2*sin(q1 - q2)**2*cos(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)**2) - m2*p1*p2*sin(q1 - q2)*sin(2*q1 - 2*q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)**2) + p1*p2*cos(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) + 2*m2*(l1**2*p2**2*(m1 + m2) - 2*l1*l2*m2*p1*p2*cos(q1 - q2) + l2**2*m2*p1**2)*sin(q1 - q2)*sin(2*q1 - 2*q2)*cos(q1 - q2)/(l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**3) - (l1**2*p2**2*(m1 + m2) - 2*l1*l2*m2*p1*p2*cos(q1 - q2) + l2**2*m2*p1**2)*cos(2*q1 - 2*q2)/(l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**2)
    d33=-p2*sin(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) + (-2*l1*l2*m2*p2*cos(q1 - q2) + 2*l2**2*m2*p1)*sin(2*q1 - 2*q2)/(2*l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**2)
    d34=-p1*sin(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) + (2*l1**2*p2*(m1 + m2) - 2*l1*l2*m2*p1*cos(q1 - q2))*sin(2*q1 - 2*q2)/(2*l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**2)
    
    d41=-2*m2*p1*p2*sin(q1 - q2)**2*cos(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)**2) - m2*p1*p2*sin(q1 - q2)*sin(2*q1 - 2*q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)**2) + p1*p2*cos(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) + 2*m2*(l1**2*p2**2*(m1 + m2) - 2*l1*l2*m2*p1*p2*cos(q1 - q2) + l2**2*m2*p1**2)*sin(q1 - q2)*sin(2*q1 - 2*q2)*cos(q1 - q2)/(l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**3) - (l1**2*p2**2*(m1 + m2) - 2*l1*l2*m2*p1*p2*cos(q1 - q2) + l2**2*m2*p1**2)*cos(2*q1 - 2*q2)/(l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**2)
    d42=-g*l2*m2*cos(q2) + 2*m2*p1*p2*sin(q1 - q2)**2*cos(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)**2) + m2*p1*p2*sin(q1 - q2)*sin(2*q1 - 2*q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)**2) - p1*p2*cos(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) - 2*m2*(l1**2*p2**2*(m1 + m2) - 2*l1*l2*m2*p1*p2*cos(q1 - q2) + l2**2*m2*p1**2)*sin(q1 - q2)*sin(2*q1 - 2*q2)*cos(q1 - q2)/(l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**3) + (l1**2*p2**2*(m1 + m2) - 2*l1*l2*m2*p1*p2*cos(q1 - q2) + l2**2*m2*p1**2)*cos(2*q1 - 2*q2)/(l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**2)
    d43=p2*sin(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) - (-2*l1*l2*m2*p2*cos(q1 - q2) + 2*l2**2*m2*p1)*sin(2*q1 - 2*q2)/(2*l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**2)
    d44=p1*sin(q1 - q2)/(l1*l2*(m1 + m2*sin(q1 - q2)**2)) - (2*l1**2*p2*(m1 + m2) - 2*l1*l2*m2*p1*cos(q1 - q2))*sin(2*q1 - 2*q2)/(2*l1**2*l2**2*(m1 + m2*sin(q1 - q2)**2)**2)
  
    """
    der_matrix stores all the derivatives in a matrix
    """
    der_matrix=np.array([[d11,d12,d13,d14],[d21,d22,d23,d24],[d31,d32,d33,d34],[d41,d42,d43,d44]])
    i=np.identity(4)
    
    res_der=i-der_matrix*h
    
    res_der_inv=np.linalg.inv(res_der)
    
    return res_der_inv
 
    """
    this function calculates the energy at every time step
    """
    
def energy_calculation(q1,q2,p1,p2,m1,m2,l1,l2,x1,x2,x3,x4):
    
    g=9.8
    P=-(m1+m2)*g*l1*cos(q1)-m2*g*l2*cos(q2)
    T=0.5*m1*l1**2*x1**2+0.5*m2*(l1**2*x1**2+l2**2*x2**2+2*l1*l2*x1*x2*cos(q1-q2))
    
    return P,T
 
def acceleration(q1,q2,m1,m2,l1,l2,x1,x2):
    g=9.81
    A=np.array([[(m1+m2)*l1,m2*l2*cos(q1-q2)],[l1*cos(q1-q2),l2]])
    B=np.array([[-(m2*l2*x2**2*sin(q1-q2)+(m1+m2)*g*sin(q1))],[-g*sin(q2)+l1*x1**2*sin(q1-q2)]])
    X = np.linalg.solve(A,B)
    
    return X[0,0],X[1,0]
    
    



   
