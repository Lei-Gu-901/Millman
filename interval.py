import numpy as np
import matplotlib.pyplot as plt
import struct

with open("stimes.dat", mode='rb') as file: # b is important -> binary
    spike = file.read()
file.close()

spike_list = []
for i in struct.iter_unpack('d', spike):
    spike_list.append(i[0])

print(len(spike_list))
# start = spike_list[10000]
# end = spike_list[90000]
# intvel = (end-start)/80000
intvel = spike_list[-1]/len(spike_list)
print(intvel)

y = np.ones(10000)
# y  = np.ones(len(spike_list))

plt.figure(figsize=(10, 8))
plt.scatter(spike_list[110000:120000], y, s=10, marker='o')
# plt.scatter(spike_list, y, s=10, marker='o')
# Show the boundary between the regions:
plt.show()
# last = spike_list[-1]
# print(last)
# interval = 0.1
# length = len(spike_list)
# # print(spike_list)
# iter_list = []
# leng = len(spike_list)
#
# nint = int(np.ceil(1000000.0/interval))
# raster = np.zeros(nint+1)
# for niter in range(length):
#     ind = int(np.floor(spike_list[niter]/interval))
#     raster[ind] += 1
#
# print(raster[:100])
# num_list = []
# counter = 0
# for niter in range(nint+1):
#     num = raster[niter]
#     if num==0 and counter>0:
#         num_list.append(int(counter))
#         counter = 0
#     counter += num
#
# print("-------------")
# print(num_list)
#
# cut = 3000
# size = [0]*cut
# for item in num_list:
#     if item <= cut:
#         size[item-1] +=1
# print(size)
#
#
# # for iter in range(1,len(size)-2):
# #     if size[iter-1]+size[iter]+size[iter+1] < 20:
# #         ind = iter
# #         break
# #
# # size =size[:ind]
# #
# # print(ind)
#
# size = np.array(size, dtype=np.float)
# norm = np.sum(size)
# size = size/norm
# xlist = []
# ylist = []
# for iter in range(len(size)):
#     xlist.append(iter + 1)
#     if size[iter]>0.0:
#         ylist.append(size[iter])
#     else: ylist.append(np.NaN)
#
# N = len(ylist)
# summation = 0.0
# lower = min(ylist)
# for iter in range(N):
#     summation += np.log(ylist[iter]/lower)
#
# print(N/summation)
#
# xlist = np.array(xlist)
# ylist = np.array(ylist)
# plt.figure(figsize=(10, 8))
# plt.plot(np.log(xlist), np.log(ylist), 'b')
# plt.show()
#
# data=np.array([xlist,ylist])