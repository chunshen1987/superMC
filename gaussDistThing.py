import numpy as np 
import matplotlib.pyplot as plt 

data = np.loadtxt('data/quarks.data')

p1 = []
p2 = []
p3 = []
j = 0
for i in range(len(data)):
	if j == 0:
		p1.append([data[i,0],data[i,1]])
	if j == 1:
		p2.append([data[i,0],data[i,1]])
	if j == 2:
		p3.append([data[i,0],data[i,1]])

	j = (j+1)%3

print (2*np.var([r[0] for r in p1])+np.var([r[0] for r in p3]))/3.0