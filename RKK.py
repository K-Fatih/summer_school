def runge_kutta(y_initial,parameters):
    import numpy as np
    from copy import copy
    from func import function_evaluation
    from func import energy_calculation,acceleration

    
    
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
    y0=np.array(y_initial).reshape(-1,1);
    
    #setting parameters of pendulum
    m1=parameters[0];
    m2=parameters[1];
    l1=parameters[2];
    l2=parameters[3];
    
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
        acc1,acc2=acceleration(q1,q2,m1,m2,l1,l2,x1,x2)
        acc_vec[0,i]=acc1
        acc_vec[1,i]=acc2
        
        i+=1
    
    total_energy=pot_ener_vec+kin_ener_vec 
    return total_energy, kin_ener_vec, pot_ener_vec,vel_vec, acc_vec,yvec,tvec
    
