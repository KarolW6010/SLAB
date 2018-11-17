import numpy as np
import matplotlib.pyplot as plt

N = 10000	#Samples
fs = 1e6
fmin = 1e6
fmax = 1.5e6
T = 1/fs

def prepare(x):
	shp = len(x)
	vecs = []

	for i in range(shp):
		testvec = "test_data[" +str(i)+"] <= 8'b" + str(x[i]) + ";"
		vecs.append(testvec)

	return vecs

k = (fmax-fmin)/(40*np.pi*T)

t = np.linspace(0,40*np.pi/fmin,N)
x = np.round(np.cos(2*np.pi*(fmin*t + k*(t**2)/2))*127)
x = np.int_(x)

plt.figure()
plt.plot(t,x)
plt.show()


binVals = []

for i in range(N):
	if x[i] > 0:
		binVals.append("0" + str(np.binary_repr(x[i], width = 7)))
	else:
		binVals.append(str(np.binary_repr(x[i], width = 8)))

outputs = prepare(binVals)

file = open("testVecs.txt","w")
for i in outputs:
	temp = i + '\n'
	file.write(temp)
file.close()