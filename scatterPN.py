#!/usr/bin/env python

#############################################################


import random
import os
import numpy
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
import matplotlib.backends.backend_pdf
from copy import copy




##############################################################
###### Function definitons
##############################################################

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap
	
##############################################################

plot_dir = "/home/ITER/oosterw/Archive/Results/hel_helwo_ITERhelcur_mis_test__2017_12_12_15_41_54"

try:
	os.chdir(plot_dir) #change this!
except Exception, e: print(e)
	#sys.exit(1)

with  open("plotdata.txt") as datafile:
	content =datafile.readlines()

	
x=[]
y=[]
z=[]
b=[]
n=[]
for i in range(len(content)):
	columns = content[i].split("\t")
	x.append(columns[0].strip())
	y.append(columns[1].strip())
	z.append(float(columns[2].strip()))	
	b.append(float(columns[3].strip())*100)
	n.append(float(columns[4].strip()))

Ntor = sorted(set(n))	

	
# #Test with random values
# for i in range(10):	
	# for j in range(10):
		# x.append(i)
		# y.append(j)
		# z.append(i*j/100.0*0.8)	
		# b.append(j)
		# n.append(random.choice([2,5,10,20]))


# #Scatter with positive and negative values	
fig = plt.figure()
plot = plt.scatter(b,y,c=numpy.sign(z),marker='o')
plt.colorbar()
plt.xlabel("Normalized beta")
plt.ylabel("Beta poloidal pedestal")
plt.title("Stable and unstable clearly distinguished")
plt.grid(b=True, which="major", linestyle="-")
plt.grid(b=True, which="minor", linestyle="--")
fig.savefig("PositiveNegativeStability.ps")


###############################################################################
###############################################################################
# Plot with different toroidal mode numbers
###############################################################################
###############################################################################

####################
#Plot with different colors
###################

N = len(Ntor)+1 
fig, ax = plt.subplots(1,1, figsize=(6,6))
tag = []
for k in range(len(x)):
	if z[k] < 0:
		tag.append(0)
	else: 
		for nij in range(len(Ntor)):
			if n[k] == Ntor[nij]:
				tag.append(nij+1) #(nij+1)/len(Ntor) )
				break
			if nij == len(Ntor)-1:
				print("Something went wrong")
				
cmap = truncate_colormap(plt.cm.gist_ncar,0.0,0.9) 
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = cmap.from_list('Custom', cmaplist, cmap.N)

bounds = numpy.linspace(0,N,N+1)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
ticklabels = ["Stable"]
for j in range(len(Ntor)):
	ticklabels.append("n="+str(Ntor[j]))
#ticklabels = ["Stable","n="+str(Ntor[0]),"n="+str(Ntor[1]),"n="+str(Ntor[2]),"n="+str(Ntor[3]),"n="+str(Ntor[4]),""]

scatplot = ax.scatter(b,y,c=tag,cmap=cmap,norm=norm)
cb = plt.colorbar(scatplot, spacing='proportional')
cb.ax.set_yticklabels('')
for j, lab in enumerate(ticklabels):
    cb.ax.text(1.2, (2 * j + 1) / (N*2.0), lab, ha='left', va='center')
cb.ax.get_yaxis().labelpad = 15

# Adding labels and grid
plt.xlabel("Normalized beta")
plt.ylabel("Beta poloidal pedestal")
plt.title("Stability comparing toroidal modes")
plt.grid(b=True, which="major", linestyle="--")
plt.grid(b=True, which="minor", linestyle=":")
ax.set_axisbelow(True)
fig.savefig("ntorStability.png")

#######################################
#### Plot with stable point as empty circles
#######################################\

xpos=[]
ypos=[]
zpos=[]
bpos=[]
npos=[]
xneg=[]
yneg=[]
zneg=[]
bneg=[]
nneg=[]
for j in range(len(z)):
	if z[j] < 0:
		xneg.append(x[j])
		yneg.append(y[j])
		zneg.append(z[j])
		bneg.append(float(b[j]))
		nneg.append(n[j])
	else:
		xpos.append(x[j])
		ypos.append(y[j])
		zpos.append(z[j])
		bpos.append(float(b[j]))
		npos.append(n[j])
		
N = len(Ntor) #Changed to get 'Stable' out of colorbar
Ncmap = N-1
fig, ax = plt.subplots(1,1, figsize=(6,6))
tag2 = []
for k in range(len(xpos)):
		for nij in range(len(Ntor)):
			if npos[k] == Ntor[nij]:
				tag2.append(nij) #(nij+1)/len(Ntor) )
				break
			if nij == len(Ntor)-1:
				print("Something went wrong")					
		
cmap = truncate_colormap(plt.cm.gist_ncar,0.0,0.9) 
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = cmap.from_list('Custom', cmaplist, cmap.N)

bounds = numpy.linspace(0,N,N+1)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

ticklabels = [] #Changed to get 'Stable' out of colorbar
for j in range(len(Ntor)):
	ticklabels.append("n="+str(Ntor[j]))

plotneg = plt.scatter(bneg,yneg,facecolors='none',edgecolors='black')
plotpos = plt.scatter(bpos,ypos,c=tag2,cmap=cmap,norm=norm)
cb = plt.colorbar(spacing='proportional')
cb.ax.set_yticklabels('')
for j, lab in enumerate(ticklabels):
    cb.ax.text(1.2, (2 * j + 1) / (N*2.0), lab, ha='left', va='center')
cb.ax.get_yaxis().labelpad = 15
plt.rc('axes',axisbelow=True)

plt.xlabel("Normalized beta")
plt.ylabel("Beta poloidal pedestal")
plt.title("Stability comparing toroidal modes")
plt.grid(b=True, which="major", linestyle="-")
plt.grid(b=True, which="minor", linestyle="--")
fig.savefig("ntorStabilityStableempty.png")



####################
#Plot with one color
###################



xpos=[]
ypos=[]
zpos=[]
bpos=[]
npos=[]
xneg=[]
yneg=[]
zneg=[]
bneg=[]
nneg=[]
for j in range(len(z)):
	if z[j] < 0:
		xneg.append(x[j])
		yneg.append(y[j])
		zneg.append(z[j])
		bneg.append(float(b[j]))
		nneg.append(n[j])
	else:
		xpos.append(x[j])
		ypos.append(y[j])
		zpos.append(z[j])
		bpos.append(float(b[j]))
		npos.append(n[j])
		
N = len(Ntor)+1 
fig, ax = plt.subplots(1,1, figsize=(6,6))
tag2 = []
for k in range(len(xpos)):
		for nij in range(len(Ntor)):
			if npos[k] == Ntor[nij]:
				tag2.append(nij+1) #(nij+1)/len(Ntor) )
				break
			if nij == len(Ntor)-1:
				print("Something went wrong")				
		
cmap = copy(plt.cm.Reds) #plt.cm.winter
palette = truncate_colormap(cmap,0.3,1.0)

plotneg = plt.scatter(bneg,yneg,c='blue',norm=norm)
plotpos = plt.scatter(bpos,ypos,c=tag2,cmap=palette,norm=norm)
cb = plt.colorbar(spacing='proportional')
cb.ax.set_yticklabels('')
for j, lab in enumerate(ticklabels):
    cb.ax.text(1.2, (2 * j + 1) / (N*2.0), lab, ha='left', va='center')
cb.ax.get_yaxis().labelpad = 15

plt.xlabel("Normalized beta")
plt.ylabel("Beta poloidal pedestal")
plt.title("Stability comparing toroidal modes")
plt.grid(b=True, which="major", linestyle="-")
plt.grid(b=True, which="minor", linestyle="--")
fig.savefig("ntorStabilityOnecolor.png")

###################################################
### Plot stable as color, unstable as contour
###################################################

fig = plt.figure()

print len(b),len(y),len(n)
#bUnique = list(set(b))
yUnique = list(set(y))
leny = len(yUnique)
lenb = int(math.ceil(float(len(b))/leny))
if (leny*lenb)>len(b):
	dif = leny*lenb-len(b)
	for k in range(dif):
		b.append(b[-1])
		n.append(n[-1])
		y.append(y[-1])
print leny, lenb
b2d = numpy.reshape(b,(lenb,leny))
y2d = numpy.reshape(y,(lenb,leny))
n2d = numpy.reshape(n,(lenb,leny))
print numpy.shape(b2d), numpy.shape(y2d), numpy.shape(n2d)
# bsteps = len(bUnique) 
# ysteps = len(yUnique)
# print len(b),len(bUnique)
# print yUnique
# print len(n),bsteps,ysteps
# n2d = [[n[ysteps*j+k] for j in range(bsteps)] for k in range(ysteps)]
# print bUnique,  yUnique, n2d

plotneg = plt.scatter(bneg,yneg,c='blue')
print Ntor
plotpos = plt.contour(b2d,y2d,n2d)

plt.clabel(plotpos, inline = 1)
plt.xlabel("Normalized beta")
plt.ylabel("Beta poloidal pedestal")
plt.title("Stability comparing toroidal modes")
plt.grid(b=True, which="major", linestyle="-")
plt.grid(b=True, which="minor", linestyle="--")
fig.savefig("ntorStabilityContourplot.png")








# Contourplot

# xUnique = list(set(x))
# yUnique = list(set(y))
# xsteps = len(xUnique) #len(y)/len(yUnique)
# ysteps = len(yUnique)
# print x,y
# print xUnique, yUnique
# print xsteps, ysteps
# print len(z)
# z2d = [[z[ysteps*j+k] for j in range(xsteps)] for k in range(ysteps)]
# print z2d

# levels = [0, 2e-1,4e-1,6e-1,8e-1] #Can be used to define contour values
# contourplot = ax.contour(xUnique,yUnique,z2d,colors='k')
# plt.clabel(contourplot, inline = 1)




###############################################################################
###############################################################################
# Plot with different colors for different growth rates
###############################################################################
###############################################################################

fig = plt.figure()
xpos=[]
ypos=[]
zpos=[]
bpos=[]
npos=[]
xneg=[]
yneg=[]
zneg=[]
bneg=[]
nneg=[]
for j in range(len(z)):
	if z[j] < 0:
		xneg.append(x[j])
		yneg.append(y[j])
		zneg.append(z[j])
		bneg.append(float(b[j]))
		nneg.append(n[j])
	else:
		xpos.append(x[j])
		ypos.append(y[j])
		zpos.append(z[j])
		bpos.append(float(b[j]))
		npos.append(n[j])
cmap = copy(plt.cm.Reds) #plt.cm.winter
palette = truncate_colormap(cmap,0.3,1.0)
#palette.set_under('0.0')

plotneg = plt.scatter(bneg,yneg,c='blue')
plotpos = plt.scatter(bpos,ypos,c=numpy.sqrt(zpos),cmap=palette)
plt.colorbar()
plt.xlabel("Normalized beta")
plt.ylabel("Beta poloidal pedestal")
plt.title("Growth rates")
plt.grid(b=True, which="major", linestyle="-")
plt.grid(b=True, which="minor", linestyle="--")
fig.savefig("growthrateStability.png")

###############################################################################
###############################################################################
# Plot growth rates for every N seperated
###############################################################################
###############################################################################

try:
	os.chdir(plot_dir) #change this!
except Exception, e: print(e)
	#sys.exit(1)

with  open("fulldata.txt") as datafile:
	content =datafile.readlines()
	
x=[]
y=[]
z=[]
b=[]
n=[]
for i in range(6,len(content)):
	columns = content[i].split("\t")
	if len(columns) == 10:
		x.append(columns[3].strip())
		y.append(columns[2].strip())
		z.append(float(columns[7].strip().split(",")[0].strip()))	
		b.append(float(columns[4].strip())*100)
		n.append(float(columns[5].strip()))
	#else:
		#print("Line structure not right for line "+str(i)+":"+columns[-1])
		
Ntor = sorted(set(n))	

newmap = "PlotsDifferentNtor"
if os.path.exists(plot_dir+"/"+newmap):
	print(plot_dir+"/"+newmap+" already exists.")
else:
	try:
		os.mkdir(plot_dir+"/"+newmap)
	except OSError:
		print("Cannot make directory"+plot_dir+"/"+newmap+" although it did not exist")
		sys.exit(1)
os.chdir(newmap)

cmap = copy(plt.cm.Reds) 
palette = truncate_colormap(cmap,0.3,1.0)

xpos=[]
ypos=[]
zpos=[]
bpos=[]
npos=[]
xneg=[]
yneg=[]
zneg=[]
bneg=[]
nneg=[]
for i in range(len(Ntor)):
	xneg.append([])
	yneg.append([])
	zneg.append([])
	bneg.append([])
	nneg.append([])
	xpos.append([])
	ypos.append([])
	zpos.append([])
	bpos.append([])
	npos.append([])
	for j in range(len(z)):
		if n[j] == Ntor[i]:
			if z[j] < 0:
				xneg[i].append(x[j])
				yneg[i].append(y[j])
				zneg[i].append(z[j])
				bneg[i].append(float(b[j]))
				nneg[i].append(n[j])
			else:
				xpos[i].append(x[j])
				ypos[i].append(y[j])
				zpos[i].append(z[j])
				bpos[i].append(float(b[j]))
				npos[i].append(n[j])
	
	fig = plt.figure()
	plotneg = plt.scatter(bneg[i],yneg[i],c='blue')
	plotpos = plt.scatter(bpos[i],ypos[i],c=numpy.sqrt(zpos[i]),cmap=palette)
	plt.colorbar()
	plt.xlabel("Normalized beta")
	plt.ylabel("Beta poloidal pedestal")
	plt.title("Growth rates n="+str(Ntor[i]))
	plt.grid(b=True, which="major", linestyle="-")
	plt.grid(b=True, which="minor", linestyle="--")
	fig.savefig("growthrateStability"+str(Ntor[i])+".png")
	plt.close(fig)
	
### Careful with antoher plot: current map and loaded file