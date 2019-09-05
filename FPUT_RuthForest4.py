from scipy import *
import matplotlib.pyplot as plt


def force(x_prv,x,x_nxt,m,n,p):
   f=m*((x_nxt-x)-(x-x_prv))+n*((x_nxt-x)**2-(x-x_prv)**2)+p*((x_nxt-x)**3-(x-x_prv)**3)
   return f


if __name__=='__main__':
   
   dt=0.125
   nmodes=10
   ak=0.0#weight for fourier mode
   Niter=433482
   N=71
   a=1.0 #spring constant of harmonic force
   b=0.0 #spring constant of phi-3 force
   c=1000 #spring constant of phi-4 force
   E=0.0

   #instantiate k1,k2,k3,k4 these are the RK4 "ghost" displacements
   
   geom=(1/(2-2**(1/3)))

   see=[0.5*geom,0.5*(1-geom),0.5*(1-geom),0.5*geom]

   dee=[geom,1-2*geom,geom,0]

   v1=zeros(N,dtype=float)
   v2=zeros(N,dtype=float)
   v3=zeros(N,dtype=float)
   v4=zeros(N,dtype=float)


  #placeholder and actual string displacements

   disp=zeros(N,dtype=float)
   disp_ghost=disp

   ghost_1=zeros(N,dtype=float)
   ghost_2=zeros(N,dtype=float)
   ghost_3=zeros(N,dtype=float)


   #position index
   pos_ind=linspace(0,N,N)

   #list of modes
   modes=linspace(1,nmodes,nmodes)
  
   #alternative just to track first given Fourier  modes
   fr=10**2                            #proportaion to time lapse per frame 
   Ntrak=Niter//fr                    #number of frames
   ttrack=linspace(0,Ntrak,Ntrak)    #time index
   foorier_sim=zeros((nmodes,Ntrak),dtype=float)      
   trek=0
   
   #running scalar for mode distribution
   ak=0.0
   ak_dot=0.0

   #initialize the displacement
   
   for i in range(N):
      disp[i]=sin((pi*i)/(N-1)) #This is just pure sine mode
      disp_ghost[i]=disp[i]
   print(disp_ghost)
   for i in range(1,N-1):
      E+=0.5*(v3[i])**2+0.5*a*((disp[i+1]-disp[i])**2+(disp[i-1]-disp[i])**2)+b*(1/3)*((disp[i]-disp[i-1])**3+(disp[i+1]-disp[i])**3)+c*0.25*((disp[i]-disp[i-1])**4+(disp[i+1]-disp[i])**4)
   print(E)

   for k in range(Niter):
      if k == 0:
         for i in range(1,N-1):
            v3[i]=force(disp[i-1],disp[i],disp[i+1],a,b,c)*(dt/2)
      for i in range(1,N-1):
         disp[i]+=dt*v3[i]*see[k%4]
      for i in range(1,N-1):
         v3[i]+=dee[k%4]*dt*force(disp[i-1],disp[i],disp[i+1],a,b,c)
         

         
#      disp_ghost=disp

      if k%fr==0 and k>0:
         for freq in range(1,nmodes):
            for i in range(N):
               ak+=disp[i]*sin(i*freq*pi/N)
               ak_dot+=v3[i]*sin(i*freq*pi/N)
            foorier_sim[freq][trek]=0.5*ak_dot**2+2*ak**2*sin(freq*pi/(2*N))**2#+(8/3)*b*ak**3*sin(freq*pi/(2*N))**3+4*c*ak**4*sin(freq*pi/(2*N))**4
            ak=0.0
            ak_dot=0.0
         trek+=1
   
   E=0.0
   for i in range(1,N-1): 
      E+=0.5*v3[i]**2+a*0.5*((disp[i+1]-disp[i])**2+(disp[i-1]-disp[i])**2)+b*(1/3)*((-disp[i-1]+disp[i])**3+(disp[i+1]-disp[i])**3)+c*0.25*((-disp[i-1]+disp[i])**4+(disp[i+1]-disp[i])**4)
   print(E)
   plt.figure()
   plt.subplot(2,1,1)
   plt.plot(pos_ind,disp)
   plt.subplot(2,1,2)
   for freq in range(nmodes):
      plt.plot(ttrack,foorier_sim[freq][:])
   plt.show()
