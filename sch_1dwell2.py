#Solve Schroedinger equation for 1D inifinte Square-Well potential. 
import matplotlib.pyplot as plt
import numpy as np
pi=2.*np.arcsin(1)
# length of the well=L
L=1.
#unit of energy EL=ground state energy of 1D infinite square well
EL=pi*pi/(2*L*L)
#integration step dx
dx=0.0001
#number of points along x from 0 to L
N=int(L/dx)
#x coordinate array and the wavefunction array
xarr=np.linspace(0,L,num=N)
psi_correct=np.zeros(N)
#input energy range(in terms of EL)
emin=0.
emax=105.
ne=1000
Earr=np.linspace(emin,emax,num=ne)
psiend=np.zeros(ne)
#define a function to compute the wavefunction at x=L using Finite-difference
def psi_at_L(energy):
    psi=np.zeros(N)
    psi[0]=0.
    psi[1]=1.
    tot=0.
    for i in range(2,N):
        psi[i]=-(1.+2.*energy*EL*dx*dx)*psi[i-2]+2.*psi[i-1]
        tot=tot+psi[i]*psi[i]*dx
    psiL=psi[N-1] 
    return psiL,tot,psi
#array to store eigenenergies ( a maximum number of iemax)
iemax=int(np.sqrt(emax))
eigenE=np.zeros(iemax)
ieg_count=0
psi_correct_all=np.zeros((iemax,N))
for k in range(0,ne):
    psiend[k],norm_const,psi_correct=psi_at_L(Earr[k])
#capture when the psiend[k] passes through zero (i.e.,find eigen energies)
    if k>0:
        a=psiend[k]
        b=psiend0
        if a<0 and b>0 or a>0 and b<0:
            ieg_count=ieg_count+1  
            eigenE[ieg_count-1]=Earr[k]
#store the correct wavefunction to print
            for j in range(0,N):
                psi_correct_all[ieg_count-1,j]=psi_correct[j]/norm_const
    psiend0=psiend[k]
print('Energy levels found:',ieg_count) 
fig1=plt.figure(1)
plt.plot(Earr,psiend)
plt.xlabel('E/E_L')
plt.ylabel('Psi(x=L)')
plt.grid()
plt.scatter(eigenE,np.zeros(ieg_count),c='r')
plt.savefig('1dsq_roots.pdf',bbox_inches = "tight")

fig2=plt.figure(2)
for j in range(0,6):   
    plt.plot(xarr,psi_correct_all[j],label='n=%s'%(j+1))
plt.xlabel('x (a.u)')
plt.ylabel('Psi(x)')
plt.grid()
plt.legend()
plt.savefig('1dsq_psi.pdf',bbox_inches = "tight")

