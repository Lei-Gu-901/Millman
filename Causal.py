import numpy as np
import matplotlib.pyplot as plt
import struct

'''
size distribution
'''
with open("num_ava.dat", mode='rb') as file: # b is important -> binary
    num = file.read()
file.close()

num_list = []
for i in struct.iter_unpack('i', num):
    num_list.append(i[0])

maxium = max(num_list)

cut = 1000
size = [0]*cut
for item in num_list:
    if item <= cut:
        size[item-1] +=1

size = np.array(size, dtype=np.float)
norm = np.sum(size)
size = size/norm
xlist = []
ylist = []
for iter in range(len(size)):
    xlist.append(iter + 1)
    if size[iter]>0.0:
        ylist.append(size[iter])
    else: ylist.append(np.NaN)

N = len(ylist)
summation = 0.0
lower = min(ylist)
for iter in range(N):
    summation += np.log(ylist[iter]/lower)

xlist = np.array(xlist)
ylist = np.array(ylist)

plt.figure(figsize=(10, 8))
plt.scatter(np.log(xlist), np.log(ylist),s=10, marker='o')
plt.show()

'''
duration distribution
'''
with open("ava_start.dat", mode='rb') as file: # b is important -> binary
    start = file.read()
file.close()

start_list = []
for i in struct.iter_unpack('d', start):
    start_list.append(i[0])

with open("ava_end.dat", mode='rb') as file: # b is important -> binary
    end = file.read()
file.close()

end_list = []
for i in struct.iter_unpack('d', end):
    end_list.append(i[0])

eint = 1.23
print(end_list[:100])
print(start_list[:100])
leng = len(end_list)

cut = 1000
size = [0]*cut
for iter in range(leng):
    item = int(np.floor((end_list[iter]-start_list[iter])/eint)) + 1
    if item <= cut:
        size[item-1] +=1
print(size)

size = np.array(size, dtype=np.float)
norm = np.sum(size)
size = size/norm
xlist = []
ylist = []
for iter in range(len(size)):
    xlist.append(iter + 1)
    if size[iter]>0.0:
        ylist.append(size[iter])
    else: ylist.append(np.NaN)


xlist = np.array(xlist)
ylist = np.array(ylist)
plt.figure(figsize=(10, 8))
plt.scatter(np.log(xlist), np.log(ylist),s=10, marker='o')
plt.show()

