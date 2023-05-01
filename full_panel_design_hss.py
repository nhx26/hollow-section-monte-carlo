import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

sigma = 247 #MPa # WOVEN ROVING ONLY
E = 12 #GPa
Ex = E
Ey = E
fs = 2.5
poisson = 0.3
k_iso = 3.5
def optimize_section(load, member_length):

	k = 1 # both ends pinned
	lc = k*member_length #effective buckling length
	Ec = 0.13
	Gc = 0.04
	poisson = 0.3
	k_iso = 3.5
	density = 1500 #kg/m3

	samplesize = 10000000

	A_req = ((fs*np.abs(load))/(sigma*1000))*10**6
	I_req = ((fs*np.abs(load)*lc**2)/((np.pi)**2*(E*10**6)))*10**12

	print("Areq = "+str(A_req) + "mm2")
	print("Ireq = "+str(I_req)+"mm4")

	hc_min = 25
	hc_max = 80

	bc_min = 25
	bc_max = 80

	t_min= 5 #esoteric interpretation of design code min 9mm split into 5 and 5
	t_max= 20

	hc_range = [30,35,40,45,50,55,60,65,70,75]#np.linspace(hc_min,hc_max,hc_max-hc_min+1)

	bc_range = [50,55,60,65,70,75]#np.linspace(bc_min,bc_max,bc_max-bc_min+1)

	t_range = np.linspace(t_min,t_max,t_max-t_min+1)

	hc = np.random.choice( hc_range, size = samplesize, replace = True, p = None)
	bc = np.random.choice( bc_range, size = samplesize, replace = True, p = None)
	t = np.random.choice( t_range, size = samplesize, replace = True, p = None)


	h= hc+2*t
	b= bc+2*t

	A = (b*h)-(bc*hc)
	Iy = ((1/12)*b**3*h)-((1/12)*bc**3*hc)
	Ix = ((1/12)*h**3*b)-((1/12)*hc**3*bc)

	sigma_actual = (load/A)*1000 #MPa

	#local stability

	sig_local = k_iso*((np.pi**2*np.sqrt(Ex*Ey))/(12*(1-poisson**2)))*((t/b)**2)
	#sig_local= 0.78*(E*Ec*Gc)**0.33333333333333333
	sig_local = sig_local*1000 #GPa to MPa
	#print("sig_local={0}".format(sig_local))
	A_delta = A - A_req
	Ix_delta = Ix - I_req
	Iy_delta = Iy - I_req


	local_delta = sig_local - (fs*sigma_actual)

	pass_index = np.where(A_delta>0)
	pass_index = pass_index[0]

	if load <0:

		pass_index = np.intersect1d(pass_index,np.where(Ix_delta>0))
		pass_index = np.intersect1d(pass_index,np.where(Iy_delta>0))
		pass_index = np.intersect1d(pass_index,np.where(local_delta>0))


	loss = A_delta[pass_index]

	optimal_index = pass_index[np.argmin(loss)]


	hc_opt  = hc[optimal_index]
	bc_opt  = bc[optimal_index]
	t_opt = t[optimal_index]
	A_opt  = A[optimal_index]


	Ix_opt = Ix[optimal_index]
	Iy_opt = Iy[optimal_index]

	member_weight = (A_opt/1000**2)*member_length*density

	return hc_opt,bc_opt,t_opt,member_weight




# Define the coordinates of the nodes
x = [0, 0, 1, 2,2]
y = [0, 1, 1, 1,0]

# Define the connections between the nodes
connections = [(0, 1), (0, 4), (0, 2), (1, 2), (2, 3),(2,4),(3,4)]

loads = [-3,85,-30,-86,-86,-30,-3]

# Create a figure and axis
fig, ax = plt.subplots()

# Plot the connections
index = 0
total_weight = 0
for connection in tqdm(connections):
    node1 = connection[0]
    node2 = connection[1]
    ax.plot([x[node1], x[node2]], [y[node1], y[node2]], 'b-',linewidth=5)

    k_iso = 3.5

    node1 = connection[0]
    node2 = connection[1]
    x_avg = (x[node1] + x[node2]) / 2
    y_avg = (y[node1] + y[node2]) / 2
    dx = x[node2] - x[node1]
    dy = y[node2] - y[node1]
    length = np.sqrt(dx**2 + dy**2)
    angle = np.arctan2(dy, dx)
    rotation = angle * 180 / np.pi
    offset = 0.02 * length # adjust the offset as desired
    dx_offset = -offset * np.sin(angle)
    dy_offset = offset * np.cos(angle)
    h,b,t,weight = optimize_section(loads[index],length)
    index = index+1
    total_weight = total_weight+weight
    ax.text(x_avg + dx_offset, y_avg + dy_offset, 'hc={0},bc={1},t={2},weight={3}'.format(h,b,t,round(weight,2)), ha='center', va='center', rotation=rotation, fontsize=18)

    #sig_local = k_iso*((np.pi**2*np.sqrt(Ex*Ey))/(12*(1-poisson**2)))*((t/b)**2)
    #print("sig_local={0}".format(sig_local))


# Set the limits of the plot
ax.set_xlim(-2, 10)
ax.set_ylim(-5, 5)

# Set the title of the plot
ax.set_title('Panel V2, sigma = {0}, E = {1}, FS = {2} Total Weight = {3} kg'.format(sigma,E,fs,round(total_weight,2)),fontsize=18)
ax.axis('equal')
ax.axis('off')

# Show the plot
plt.show()
