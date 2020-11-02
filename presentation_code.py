"""
This file is the main file

for part d of the task sheet, we just have to call the function part_d with 
initial conditions and paraneters as two lists as input

for part f of the task sheet, we can call part_f_explicit, part_f_implicit
or part_f_rk method. for these three functions we can input our initial
conditions [theta1_initial, theta2_initial, p1_initial, p2_initial]
in the form of a list and also the parameters of the pendulum [m1,m2,l1,l2]
"""



"""
Solution of Part d
"""
def part_d(y_initial,parameters):
    
    import numpy as np
    from explicit_euler import explicit_euler
    from implicit_euler import implicit_euler
    from RKK import runge_kutta
    import matplotlib.pyplot as plt
    
    

    total_energy_1,kin_ener_vec_1,pot_ener_vec_1,vel_vec_1,acc_vec_1,yvec_1,tvec=explicit_euler(y_initial,parameters)
    
    dev_x=tvec
    plt.figure(0)
    plt.plot(dev_x,total_energy_1)
    plt.plot(dev_x,kin_ener_vec_1)
    plt.plot(dev_x,pot_ener_vec_1)
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.legend(["Total Energy","Kinetic Energy","Potential Energy"])
    plt.title("Conservation of Energy Over Time Using Explicit Euler")
    mean_energy_1=np.mean(total_energy_1)
    deviation_mean_1=[x-mean_energy_1 for x in total_energy_1]
    plt.figure(1)
    plt.plot(dev_x,deviation_mean_1)
    plt.xlabel("Deviation from mean value")
    plt.ylabel("Time")
    plt.title("Deviation of Total Energy from the Mean Value Over Time for Explicit Euler")


   
    total_energy_2,kin_ener_vec_2,pot_ener_vec_2,yvec_2,vel_vec_2,acc_vec_2,tvec=implicit_euler(y_initial,parameters)   
    dev_x=tvec
    plt.figure(2)
    plt.plot(dev_x,total_energy_2)
    plt.plot(dev_x,kin_ener_vec_2)
    plt.plot(dev_x,pot_ener_vec_2)
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.legend(["Total Energy","Kinetic Energy","Potential Energy"])
    plt.title("Conservation of Energy Over Time Using Implicit Euler")
    mean_energy_2=np.mean(total_energy_2)
    deviation_mean_2=[x-mean_energy_2 for x in total_energy_2]
    plt.figure(3)
    plt.plot(dev_x,deviation_mean_2)
    plt.xlabel("Deviation from mean value")
    plt.ylabel("Time")
    plt.title("Deviation of Total Energy from the Mean Value Over Time for Implicit Euler")



    
    
    total_energy_3,kin_ener_vec_3,pot_ener_vec_3,vel_vec_3,acc_vec_3,yvec_3,tvec=runge_kutta(y_initial,parameters)
    dev_x=tvec
    plt.figure(4)
    plt.plot(dev_x,total_energy_3)
    plt.plot(dev_x,kin_ener_vec_3)
    plt.plot(dev_x,pot_ener_vec_3)
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.legend(["Total Energy","Kinetic Energy","Potential Energy"])
    plt.title("Conservation of Energy Over Time Using Runge-Kutta")
    plt.title("Conservation of Energy Over Time Using Implicit Euler")
    mean_energy_3=np.mean(total_energy_3)
    deviation_mean_3=[x-mean_energy_3 for x in total_energy_3]
    plt.figure(5)
    plt.plot(dev_x,deviation_mean_3)
    plt.xlabel("Deviation from mean value")
    plt.ylabel("Time")
    plt.title("Deviation of Total Energy from the Mean Value Over Time for RK Method")
    



"""
Solution of part f

"""    
def part_f_explicit(y_initial,parameters):
    

    from explicit_euler import explicit_euler
    import matplotlib.pyplot as plt
    
    total_energy_1,kin_ener_vec_1,pot_ener_vec_1,vel_vec_1,acc_vec_1,yvec_1,tvec=explicit_euler(y_initial,parameters)
    dev_x=tvec
    plt.figure(3)
    plt.plot(dev_x,yvec_1[0,:])
    plt.plot(dev_x,yvec_1[1,:])
    plt.xlabel("Time")
    plt.ylabel("Position")
    plt.legend(["Mass 1","Mass 2"])
    plt.title("Position of Mass 1 and Mass 2 Over Time Using Explicit Method")

    plt.figure(4)
    plt.plot(dev_x,vel_vec_1[0,:])
    plt.plot(dev_x,vel_vec_1[1,:])
    plt.xlabel("Time")
    plt.ylabel("Velocity")
    plt.legend(["Mass 1","Mass2"])
    plt.title("Velocity of Mass 1 and Mass 2 Over Time Using Explicit Method")
    

    plt.figure(5)
    plt.plot(dev_x,acc_vec_1[0,:])
    plt.plot(dev_x,acc_vec_1[1,:])    
    plt.xlabel("Time")
    plt.ylabel("Acceleration")
    plt.legend(["Mass 1","Mass2"])
    plt.title("Acceleration of Mass1 and Mass 2 Over Time Using Explicit Method")
 

    
def part_f_implicit(y_initial,parameters):
    
    
    from implicit_euler import implicit_euler
    import matplotlib.pyplot as plt
    
    total_energy_1,kin_ener_vec_1,pot_ener_vec_1,vel_vec_1,acc_vec_1,yvec_1,tvec=implicit_euler(y_initial,parameters)
    dev_x=tvec
    plt.figure(3)
    plt.plot(dev_x,yvec_1[0,:])
    plt.plot(dev_x,yvec_1[1,:])
    plt.xlabel("Time")
    plt.ylabel("Position")
    plt.legend(["Mass 1","Mass 2"])
    plt.title("Position of Mass 1 and Mass 2 Over Time Using Implicit Method")

    plt.figure(4)
    plt.plot(dev_x,vel_vec_1[0,:])
    plt.plot(dev_x,vel_vec_1[1,:])
    plt.xlabel("Time")
    plt.ylabel("Velocity")
    plt.legend(["Mass 1","Mass2"])
    plt.title("Velocity of Mass 1 and Mass 2 Over Time Using Implicit Method")
    

    plt.figure(5)
    plt.plot(dev_x,acc_vec_1[0,:])
    plt.plot(dev_x,acc_vec_1[1,:])    
    plt.xlabel("Time")
    plt.ylabel("Acceleration")
    plt.legend(["Mass 1","Mass2"])
    plt.title("Acceleration of Mass 1 and Mass 2 Over Time Using Implicit Method")

def part_f_rk(y_initial,parameters):
    
    

    from RKK import runge_kutta
    import matplotlib.pyplot as plt
    
    total_energy_1,kin_ener_vec_1,pot_ener_vec_1,vel_vec_1,acc_vec_1,yvec_1,tvec=runge_kutta(y_initial,parameters)
    dev_x=tvec
    plt.figure(3)
    plt.plot(dev_x,yvec_1[0,:])
    plt.plot(dev_x,yvec_1[1,:])
    plt.xlabel("Time")
    plt.ylabel("Position")
    plt.legend(["Mass 1","Mass 2"])
    plt.title("Position of Mass 1 and Mass 2 Over Time Using RK Method")

    plt.figure(4)
    plt.plot(dev_x,vel_vec_1[0,:])
    plt.plot(dev_x,vel_vec_1[1,:])
    plt.xlabel("Time")
    plt.ylabel("Velocity")
    plt.legend(["Mass 1","Mass2"])
    plt.title("Velocity of Mass 1 and Mass 2 Over Time Using RK Method")
    

    plt.figure(5)
    plt.plot(dev_x,acc_vec_1[0,:])
    plt.plot(dev_x,acc_vec_1[1,:])    
    plt.xlabel("Time")
    plt.ylabel("Acceleration")
    plt.legend(["Mass 1","Mass2"])
    plt.title("Acceleration of Mass1 and Mass 2 Over Time Using RK Method")

part_d([1.04,1.04,0,0],[1,1,1,1])   