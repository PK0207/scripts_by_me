import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

# Some functions that are used to scale the density
def Hp(r):
  return np.sqrt(k*T(r)/(mu*mH)*r**3/GM)

def T(r):
  return T0 * (r/1.0)**(-1.0*q)

def Sigma(r):
  return Sigma_C * (r/r_c)**(-1.0*g)*np.exp(-1.0*(r/r_c)**(2-g))

def rho(r,z):
  H = Hp(r)
  S = Sigma(r)
  return S/(2*np.pi)**0.5 * np.exp(-0.5*(z/H)**2)/H


# Ranges for the plot
r0, r1 = 0.01, 10 # in AU
z0, z1 = 0,6
r = np.arange(r0, r1, 0.02) # adjust the last value for speed vs resolution
z = np.arange(z0, z1 ,0.02) # (higher values -> faster, smaller values -> nicer)

#print len(z)

# Some generic parameters (not physical)
T0 = 3500.0
Sigma_C = 1.0
r_c = 1.0
q = -0.1
g = 0.0
GM = 1.
mu = 1.
mH = 1.
k = 1.e-8
sig_0 = 1.1*10**-18


# generate image
im = np.zeros((len(r),len(z)))
for i in range(len(r)):
  for j in range(len(z)):
    rr = rho(r[i],z[j])
    im[i,j] = np.sqrt(rr*sig_0)
    if im[i,j] < 0.0:
      im[i,j] = 0.0


imag = np.transpose(np.log(im))
imag = (imag+370)/30.

# Display the "image"
plt.imshow(imag,aspect='auto', origin='lower', extent=(r0,r1,z0,z1))
plt.colorbar(spacing='proportional',aspect=20,ticks=[0,2,4,6,8,10],pad=0.01,label=r"H$_2$ optical depth ($\tau_{\lambda}$)")

plt.plot(r,34.*Hp(r), zorder=1, alpha=0.4, lw=2, ls='--', color='k')
plt.xlabel("r (AU)") # in fact, should be AU for the correct parameter choices.
plt.ylabel("z (AU)")
#plt.scatter([5],[2.1], marker='s',  color='w')
plt.annotate("", xy=(0.1,0), xytext=(5,2.6), arrowprops=dict(arrowstyle="<-", lw=3, color='purple'), color='purple')
plt.annotate(r"Ly$\alpha$", xy=(1.7,2.3), fontsize=20, color='purple', zorder=10)
plt.annotate(r"[$\nu$, $J$] $\rightarrow$ [$\nu$',$J$']", xy=(0.3,2), fontsize=18, color='purple', zorder=10)

plt.annotate("", xy=(5,2.6), xytext=(7.0,5), arrowprops=dict(arrowstyle="<-", lw=3, color='r'), color='r')
plt.annotate(r"H$_2$", xy=(3.5,4.0), fontsize=20, zorder=10, color='r')
plt.annotate(r"[$\nu$', $J$'] $\rightarrow$ [$\nu$'',$J$'']", xy=(1.8,3.7), fontsize=18, color='r', zorder=10)

plt.annotate(r"$\tau _{\lambda}'(r, z) \sim 1$", xy=(6.5,4.4), rotation=55, color='k', fontsize=20)
#plt.annotate(r"$\tau _{\lambda}'(r, z)$", xy=(3,2.4), rotation=-70, color='r', fontsize=20)

plt.plot([5,5],[2.7,2.2], lw=5, color='w')
plt.annotate("$\Delta z$", xy=(5.1,2.4), color='w', fontsize=30)

plt.annotate(r"$F_{Ly\alpha}$", xy=(0.55,1.0), color='purple', fontsize=20)
plt.annotate(r"$F_{H_2}$", xy=(6.5,5.1), color='r', fontsize=20)
plt.xlim(r0,r1)
plt.ylim(z0,z1)

plt.show()
#filename = 'C:\Users\keho8439\Documents\Comps2\models\model_RT2.ps'
#plt.savefig(filename,orientation='landscape')#,format=ps)
