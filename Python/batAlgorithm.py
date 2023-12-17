import numpy as np

def initializa(n,d,n_gen,lowers,uppers,optim):

    best_var_local=np.zeros([n_gen,d]) #initialize best local variable matrix
    best_var_global=np.zeros([n_gen,d]) #initialize best global variable matrix
    fitness=np.zeros(n) #fitness function
    solution=np.zeros([n,d]) #solutions for each variable
    obj_function=np.zeros(n) #objective function

    for bat in range(0,n):
        for prj_var in range(0,d):
            solution[bat,prj_var]=lowers[prj_var]+(uppers[prj_var]-lowers[prj_var])*np.random.rand()
                    
        #Compute value of objective function
        obj_function[bat]=solution[bat,0]**2+solution[bat,0]-1
    
    for bat in range(0,n):
        fitness[bat] = fitness[bat] + obj_function[bat]

    local_min, local_max=np.zeros(n_gen), np.zeros(n_gen)
    global_min, global_max= np.zeros(n_gen), np.zeros(n_gen)
    local_min[0], local_max[0]= fitness[0], fitness[0]

    for bat in range(1,n):
        if optim and fitness[bat]<=fitness[bat-1]:
            local_min[0]=fitness[bat]
            best_var_local[0,:]=solution[bat,:]

        if (not optim) and fitness[bat]>=fitness[bat-1]:
            local_max[0]=fitness[bat]
            best_var_local[0,:]=solution[bat,:]

    global_min[0]=local_min[0]
    global_max[0]=local_max[0]
    best_var_global[0,:]=best_var_local[0,:]

    return best_var_local, best_var_global, fitness, solution, obj_function, global_min, global_max, local_min, local_max

def bat_development(n,n_gen,d,A,r,f_peri,f_perf,hibrid,adaptive,Qmin,Qmax,epsilon,obj_function,global_min,local_min,local_max,global_max,best_var_global,alpha,gamma):
    
    Q=np.zeros(n)

    v=np.zeros([n,d])
    S=np.zeros([n,d])
    vec=np.zeros([n,d])
    A[1,:]=A[0,:]
    r[1,:]=r[0,:]
    BESTvar=np.zeros(d)

    for gen in range(1,n_gen):

        #if hibrid bat with DE
        if hibrid==1:
            f_per=f_peri+((f_perf-f_peri)/n_gen)*(gen-2)
        else:
            f_per=f_peri
    
        for i in range(0,n):
            
            Q[i]=Qmin+(Qmax-Qmin)*np.random.rand()
            
            for j in range(0,d):
           
                v[i,j]=v[i,j]+(solution[i,j]-best_var_global[i,j])*Q[i]
                S[i,j]=S[i,j]+v[i,j]
       
        if hibrid==0:
        
                for i in range(0,n):
                    raleat=np.random.rand()
                    bat=np.random.randint(n)
                    while bat==i:
                        bat=np.random.randint(n)
                    
                    if raleat>r(gen,bat):
                        
                        for a in range(0,d):
                            vec[0,a]=np.random.rand()
                        
                        S[i,:]=best_var_global[gen-1,:]+vec[0,:]*epsilon
                    
                for i in range(0,n):
                    for j in range(0,d):
                        
                        if S[i,j]>uppers[j]:
                            S[i,j]=uppers[j]
                        elif S[i,j]<lowers[j]:
                            S[i,j]=lowers[j]
                                    
                
        elif hibrid==1:
                
                trial=np.zeros([n,d])
                
                for k in range(0,n):
                    for m in range(0,d):
                        ind2=np.random.randint(n)
                        if ind2==k: np.random.randint(n)
                        ind3=np.random.randint(n)
                        if ind3==k: np.random.randint(n)
                        ind4=np.random.randint(n)
                        if ind4==k: np.random.randint(n)
                        ind5=np.random.randint(n)
                        if ind5==k: np.random.randint(n)
                    
                    if n<=6:
                        trial[k,m]=best_var_global[gen-1,m]+f_per*(S[ind2,m]-S[ind3,m])
                    else:
                        trial[k,m]=best_var_global[gen-1,m]+f_per*(S[ind2,m]-S[ind3,m]+S[ind4,m]-S[ind5,m])

                    rand=np.random.rand()
                    S[k,m]=rand*best_var_global[gen-1,m]+(1-rand)*trial[k,m];
                    
                for i in range(0,n):
                    for j in range(0,d):
                        
                        if S[i,j]>uppers[j]:
                            S[i,j]=uppers[j]
                        elif S[i,j]<lowers[j]:
                            S[i,j]=lowers[j]
    
        fitnessnew=np.zeros(n) #fitness function
        for i in range(0,n):
            
            obj_function[i]=solution[i,0]**2+solution[i,0]-1
            fitnessnew[i]=fitnessnew[i]+obj_function[i]

        #Update solution
        if optim:
            global_min[gen]=fitnessnew[0]
            
            for j in range(d):
                 best_var_global[gen,j]=S[0,j]
        
            if gen==1: 
                BEST=global_min[gen]
                for j in range(d):
                    BESTvar[j]=S[0,j]
            
            
        elif ~optim:
            best_var_global[gen,j]=fitnessnew[0]

            for j in range(d):
                best_var_global[gen,j]=S[0,j]
                 
            if gen==1:
                BEST=global_max[gen]
                for j in range(d):
                    BESTvar[j]=S[0,j]
      
      
    for i in range(1,n):
        for j in range(d):
            
            if optim and fitnessnew[i]<=global_min[gen]:
                best_var_local[gen,j]=S[i,j]
                local_min[gen]=fitnessnew[i]
                solution[i,j]=S[i,j]
                fitness[i]=fitnessnew[i]
                best_var_global[gen,j]=S[i,j]

                if local_min[gen]<BEST:
                     global_min[gen]=local_min[gen]
                     best_var_global[gen,j]=S[i,j]
                     BEST=global_min[gen]

                     for _d in range(d):
                        BESTvar[_d]=S[i,_d]
                    
                     if gen<n_gen-1:
                        A[gen+1,i]=alpha*A[gen,i]
                        r[gen+1,i]=r[gen,i]*(1-np.e**(-gamma*gen))                        
                else:
                    
                    global_min[gen]=fitnessnew[i]
                    best_var_global[gen,j]=best_var_global[gen-1,j]
                    BEST=global_min[gen]
                            
                    if gen<n_gen-1:
                        A[gen+1,:]=alpha*A[gen,:]
                        r[gen+1,:]=r[gen,:] 
             
            elif ~optim and fitnessnew[i]>=global_max[gen]:
                best_var_local[gen,j]=S[i,j]
                local_max[gen]=fitnessnew[i]
                solution[i,j]=S[i,j]
                fitness[i]=fitnessnew[i]
                best_var_global[gen,j]=S[i,j]

                if local_max[gen]>BEST:
                     global_max[gen]=local_max[gen]
                     best_var_global[gen,j]=S[i,j]
                     BEST=global_max[gen]

                     for _d in range(d):
                        BESTvar[_d]=S[i,_d]
                    
                     if gen<n_gen-1:
                        A[gen+1,i]=alpha*A[gen,i]
                        r[gen+1,i]=r[gen,i]*(1-np.e**(-gamma*gen))                        
                else:
                    
                    global_max[gen]=fitnessnew[i] #maxG(1,t)=maxG(1,t-1);??
                    best_var_global[gen,j]=best_var_global[gen-1,j]
                    BEST=global_max[gen]
                            
                    if gen<n_gen-1:
                        A[gen+1,:]=alpha*A[gen,:]
                        r[gen+1,:]=r[gen,:] 
  
# In case of having BAT hibrid DE

f_peri=1; #initial frequency 
f_perf=.3 #final frequency
adaptive=0 #adaptive approach (0: no; 1: yes)
hibrid=1 #hibrid approach (0: no; 1: yes)

# The entrance also can be interactive

n=4; #population dimension
n_gen=5 #number of generations
A0=.9 #loudness, constant or growning)
r0=.4 #pulse rate, constant ou decreasing)
alpha=.9
gama=.9

optim=True; #if we want to minimize of maximize (true: minimize; false: maximaze)

# Frequency range defines the scaling, so it might require adjustments

Qmin,Qmax= 0, 3 #min and max frqs
epsilon=1/1000; #scaling factor
d=1 #dimension of project variables

#Lower and upper limites of the project's variables' vectors

lowers=np.zeros(d)
uppers=np.zeros(d)

lowers[:]=-5
uppers[:]=10

# Arrays' Initialization

A=np.zeros([n_gen,n]) #loudness matrix
r=np.zeros([n_gen,n]) #pulse rate matrix

for i in range(0,n):
    
    A[0,i]=A0+np.random.rand()*.1
    r[0,i]=r0+np.random.rand()*.1

best_var_local, best_var_global, fitness, solution, obj_function, global_min, global_max, local_min, local_max = initializa(n,d,n_gen,lowers,uppers,optim)

print("********************************************")
print("Initialization values (first generation):")
print(f"Best local vars: {best_var_local[0,:]}")
print(f"Best global vars: {best_var_global[0,:]}")
print("--------------------------------------------")
print(f"Positions for each bat:")
for bat in range(0,n):
    print(f"Bat {bat+1} ---> {obj_function[bat]}")
print(f"Project variables for each bat:")
for bat in range(0,n):
    print(f"Bat {bat+1} ---> {solution[bat,:]}")
print("********************************************")

bat_development(n,n_gen,d,A,r,f_peri,f_perf,hibrid,adaptive,Qmin,Qmax,epsilon,obj_function,global_min,local_min,local_max,global_max,best_var_global,alpha,gama)

print("********************************************")
print("Initialization values (Last generation):")
print(f"Best local vars: {best_var_local[0,:]}")
print(f"Best global vars: {best_var_global[0,:]}")
print("--------------------------------------------")
print(f"Positions for each bat:")
for bat in range(0,n):
    print(f"Bat {bat+1} ---> {obj_function[bat]}")
print(f"Project variables for each bat:")
for bat in range(0,n):
    print(f"Bat {bat+1} ---> {solution[bat,:]}")
print("********************************************")