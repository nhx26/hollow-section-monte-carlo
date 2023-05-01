import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def plot_section(h,b,t):
	fig, ax = plt.subplots()
	ax.add_patch(Rectangle(((b/2)-t, -h/2),t,h-t,facecolor="black"))
	ax.add_patch(Rectangle((-b/2,-h/2),t,h-t,facecolor="black"))
	ax.add_patch(Rectangle((-b/2,(h/2)-t),b,t,facecolor="black"))
	plt.axis('equal')


load = float(input("Enter axial load(kN): "))#kN
member_length = float(input("Enter member length(m): "))


k = 1 # both ends pinned


fs = 2.5# factor of safety
lc = k*member_length #effective buckling length

sigma = 247 #MPa # WOVEN ROVING ONLY
E = 12 #GPa
Ex = 12
Ey = 12
Ec = 0.13
Gc = 0.04
poisson = 0.3
k_iso = 1.33
density = 1570 #kg/m3

samplesize = 10000000

A_req = ((fs*np.abs(load))/(sigma*1000))*10**6
I_req = ((fs*np.abs(load)*lc**2)/((np.pi)**2*(E*10**6)))*10**12

print("Areq = "+str(A_req) + "mm2")

if load < 0:
	print("Ireq = "+str(I_req)+"mm4")

h_min = 50
h_max = 300

b_min = 50
b_max = 300

t_min= 5
t_max= 20

h_range = np.linspace(h_min,h_max,h_max-h_min+1)

b_range = np.linspace(b_min,b_max,b_max-b_min+1)

t_range = np.linspace(t_min,t_max,t_max-t_min+1)

h = np.random.choice( h_range, size = samplesize, replace = True, p = None)
b = np.random.choice( b_range, size = samplesize, replace = True, p = None)
t = np.random.choice( t_range, size = samplesize, replace = True, p = None)


A = (b*h)-((b-2*t)*(h-t))



Iy = ((1/12)*b**3*h)-((1/12)*((b-2*t)**3)*(h-t))

Ipx = ((1/3)*b*h**3)-((1/3)*(b-2*t)*(h-t)**3)

yc = ((h**2*t)+(t*(h-(t/2))*(b-2*t)))/A

Ix = Ipx -(A*yc**2)

#local stability

sig_local = k_iso*((np.pi**2*np.sqrt(Ex*Ey))/(12*(1-poisson**2)))*((t/b)**2)
#sig_local= 0.78*(E*Ec*Gc)**0.33333333333333333
sig_local = sig_local*1000 #GPa to MPa

print(sig_local)


A_delta = A - A_req
Ix_delta = Ix - I_req
Iy_delta = Iy - I_req

local_delta = sig_local - (fs*sigma)

pass_index = np.where(A_delta>0)
pass_index = pass_index[0]

if load <0:

	pass_index = np.intersect1d(pass_index,np.where(Ix_delta>0))
	pass_index = np.intersect1d(pass_index,np.where(Iy_delta>0))
	pass_index = np.intersect1d(pass_index,np.where(local_delta>0))


loss = A_delta[pass_index]

optimal_index = pass_index[np.argmin(loss)]


h_opt  = h[optimal_index]
b_opt  = b[optimal_index]
t_opt = t[optimal_index]
A_opt  = A[optimal_index]


Ix_opt = Ix[optimal_index]
Iy_opt = Iy[optimal_index]

member_weight = (A_opt/1000**2)*member_length*density

print("h: "+str(h_opt)+"mm")
print("b: "+str(b_opt)+"mm")
print("t: "+str(t_opt)+"mm")

print("Area: "+str(A_opt)+"mm2")
print("Ix: "+str(Ix_opt)+"mm4")
print("Iy: "+str(Iy_opt)+"mm4")

print("Weight: "+str(member_weight)+"kg")

plot_section(h_opt,b_opt,t_opt)
plt.show()