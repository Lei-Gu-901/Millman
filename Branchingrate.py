import numpy as np
import matplotlib.pyplot as plt
import struct

with open("num_branch.dat", mode='rb') as file: # b is important -> binary
    branch = file.read()
file.close()

with open("num_spike.dat", mode='rb') as file: # b is important -> binary
    spike = file.read()
file.close()

# exclude the receipt and transmission layers (first six layers)
n = 6
start = n*81
end = (10+n)*81
length = 10+n

brach_list = []
for i in struct.iter_unpack('i', branch):
    brach_list.append(i[0])

spike_list = []
for i in struct.iter_unpack('i', spike):
    spike_list.append(i[0])

spike = np.array(spike_list)
branch = np.array(brach_list)

leng = len(spike_list)
# in the rare case that an isolated neuron exist so that it spike number is 0
# then set it as 1 to prevent the error of divided by 0
for iter in range(leng):
    if spike[iter] == 0:
        spike[iter] = 1

branchingrate = branch/spike
print(np.mean(branchingrate[start:end]))
