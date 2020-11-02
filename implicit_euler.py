def implicit_euler(y_initial,parameters):
    import numpy as np
    from copy import copy
    from func import function_evaluation, der_func,energy_calculation,acceleration
    
    """
    t0=starting time, T=maximum time, h=step size, y0=inital condition array
    tvec stores value of time, yvec stores yj at each time step. yj is a 4x1 array
    with the first entry theta1,second entry theta2, third entry p1
     and the fourth entry q2
     
    tol is the tolerance of the residuum for the newton rapshon scheme 
    """
    t0=0;
    T=10;
    h=0.001;
    y0=np.array(y_initial);
    
    #setting parameters of pendulum
    m1=parameters[0];
    m2=parameters[1];
    l1=parameters[2];
    l2=parameters[3];
    
    num_steps=round((T-t0)/h);
    
    #values of the four variables will be stored in yvec for each time step and
    #the corresponding time step will be stored in tvec
    
    yvec = np.zeros((4,num_steps))
    tvec=np.zeros(num_steps)
    
    
    yj=y0.copy()
    tj=copy(t0)
    
    
    
    #getting the values of RK method from the butcher table
    #m indicates the order of the RK method
    #num_eqts indicate the number of equations in the RHS function vector
    
    tol=1e-10;
    #counter refers to the current time step number
    counter=0;
    while counter<num_steps:
        yvec[:,counter]=yj
        tvec[counter]=tj
        
    
        
        """
        Implicit euler method starts from here
        """
        #initial guess of yj1=yj
        yj1_old=yj.copy()
        q1=yj1_old[0]
        q2=yj1_old[1]
        p1=yj1_old[2]
        p2=yj1_old[3]
        
        """
        function_evaluation is also defined  in func file
        to understand implicit euler scheme, please refer to the derivation that I
        uploaded
        
        """
        
            
        """
        x1,x2,x3,x4 are the four hamiltonian equations
        """
            
        x1,x2,x3,x4=function_evaluation(q1,q2,p1,p2,m1,m2,l1,l2)
        
    
    
        function_val_old=np.array([x1,x2,x3,x4])
        yj1_new=yj+function_val_old*h
    
        residum=yj1_old-yj-function_val_old*h
        
        error=abs(np.average(residum))
        iter=0
        
        """
        newton raphson scheme starts here
        """
        while error > tol:
            #function_evaluation call gives unc
            #function_derivative call gives der_inv
    
            """
            calculating the values of tangent
            
            """
            
            res_der_inv=der_func(q1,q2,p1,p2,m1,m2,l1,l2)
            
            """
            updating the value
            
            """
            
            yj1_new=yj1_old-np.dot(res_der_inv,residum)
            
            """
            calculating the new value of residum
            """
            
            yj1_old=yj1_new.copy()
            q1=yj1_old[0]
            q2=yj1_old[1]
            p1=yj1_old[2]
            p2=yj1_old[3]
            
            x1,x2,x3,x4=function_evaluation(q1,q2,p1,p2,m1,m2,l1,l2)
            function_val_old=np.array([x1,x2,x3,x4])
            residum=yj1_old-yj-function_val_old*h
            
            
            
            error=abs(np.average(residum))
            """
            if number of iterations become too large, we break outta loop here
            """
        
            if iter>20:
                print(iter)
                
                break
            
    
            iter+=1
        
        tj+=h
        yj=yj1_new.copy()
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
        acc1,acc2=acceleration(q1,q2,m1,m2,l1,l2,x1,x2)
        acc_vec[0,i]=acc1
        acc_vec[1,i]=acc2
        
        i+=1
    
    total_energy=pot_ener_vec+kin_ener_vec  
    return total_energy, kin_ener_vec, pot_ener_vec,vel_vec, acc_vec,yvec,tvec
       
