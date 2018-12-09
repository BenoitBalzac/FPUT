from scipy import *
import matplotlib.pyplot as plt


def force(x_prv,x,x_nxt,a,b):
   f=0.0
   f=a*(x_prv-x+x_nxt-x)+b*((x_prv-x)**3+(x_nxt-x)**3)
   return f


if __name__=='__main__':
   
   dt=0.1
   nmodes=10
   ak=0.0#weight for fourier mode
   Niter=328078
   N=100+1
   a=1 #spring constant of harmonic force
   b=10**3 #spring constant of anharmonic force
   
   #instantiate veloities,displacements, its ghost and the force array

   vel=zeros(N,dtype=float)
   vel_ghost=zeros(N,dtype=float)
   disp=zeros(N,dtype=float)
   disp_ghost=disp
   foor=zeros(N,dtype=float)

   pos_ind=linspace(0,N,N)

   #storage array for mode distribution
   foorier=zeros(nmodes,dtype=float)
   modes=linspace(1,nmodes,nmodes)
   #running scalar for mode distribution
   ak=0.0
   ak_dot=0.0
#initialize the displacement

   for i in range(N):
      disp[i]=sin((2*pi*i)/(N-1)) #This is just pure sine mode

   for k in range(Niter):
      for i in range(1,N-1):
         foor[i]= force(disp_ghost[i-1],disp_ghost[i],disp[i+1],a,b)
      
      disp=disp_ghost+vel_ghost*dt+0.5*foor*(dt**2)
      disp_ghost=disp
      vel=vel_ghost+foor*dt
      
   for mode in range(1,nmodes):
      for i in range(N):
         ak+=disp[i]*sin((i*mode*pi)/N)
         ak_dot+=vel[i]*sin((i*mode*pi)/N)
      foorier[mode-1]=0.5*(ak_dot**2)+2*(ak**2)*(sin((mode*pi)/(2*N))**2)
      ak=0.0
      ak_dot=0.0


plt.figure()
plt.subplot(2,1,1)
plt.plot(pos_ind,disp)
plt.subplot(2,1,2)
plt.plot(modes,foorier)
plt.show()
