import numpy as np
import matplotlib.pyplot as plt
import struct

with open("stimes.dat", mode='rb') as file: # b is important -> binary
    timing = file.read()
file.close()

with open("sind.dat", mode='rb') as file: # b is important -> binary
    index = file.read()
file.close()


'''
Spikes are delta functions in time; they are broadened with Gaussian.
Then avalanches are divided into left and right inclined ones according to the peak location.
Avalanches of size range [3,104] are analyzed.
The plot is about size range [7,100], since an average over nearby 8 size is taken for a certain size.   
'''


time_list = []
for i in struct.iter_unpack('d', timing):
    time_list.append(i[0])


ind_list = []
for i in struct.iter_unpack('i', index):
    ind_list.append(i[0])


dense = 1                     # subsampling per dense neuron; dense=1 means full sampling
width_factor = 2              # scaling factor to the original time-bin width

etime = []
leng = len(ind_list)
for iter in range(leng):
    tind = ind_list[iter]
    if tind%dense== 0:
        etime.append(time_list[iter])


eint = width_factor*0.123*dense
eleng = len(etime)
llim = 3
slim = 104

srate = np.zeros(slim-llim-7, dtype=np.float)
vx = np.linspace(-5, 105, 200)


def gaussian(mu, gamma):
    vmu = vx - mu
    norm = 1.0 / np.sqrt(2.0 * np.pi) / gamma
    rv = norm * np.exp(-0.5 * vmu * vmu / (gamma * gamma))
    return rv


elist = []
numlist = []
tem = []
etime = np.array(etime)
bin = np.floor(etime/eint)
start = bin[0]
tem.append(etime[0])

for niter in range(1,eleng):
    ptime = bin[niter]
    if ptime-start< 1.1:
        tem.append(etime[niter])
    else:
        elist.append(tem)
        numlist.append(len(tem))
        tem = []
    start = ptime


spkies = [[] for _ in range(slim-llim+1) ]
llen = len(numlist)
for iter in range(llen):
    tlen = numlist[iter]
    if llim-1<tlen<slim+1:
        spkies[tlen-llim].append(elist[iter])


his = [[],[]]

for liter in range(8):
    snl = 0
    snr = 0
    savas = spkies[liter]
    gamma = 200/(liter+llim)
    for ava in savas:
        tlen = len(ava)
        arr = np.array(ava)
        arr = arr - arr[0]
        arr = arr*100/arr[-1]
        wave = np.zeros(200, dtype=np.float)
        for iter in range(tlen):
            wave += gaussian(arr[iter], gamma)
        wave = wave/tlen
        mind = np.argmax(wave)
        if mind < 100:
            snl += 1
        else:
            snr += 1
    his[0].append(snl)
    his[1].append(snr)

for liter in range(8, slim-llim+1):

    snl = 0
    snr = 0
    savas = spkies[liter]
    gamma = 200/(liter+llim)
    for ava in savas:
        tlen = len(ava)
        arr = np.array(ava)
        arr = arr - arr[0]
        arr = arr*100/arr[-1]
        wave = np.zeros(200, dtype=np.float)
        for iter in range(tlen):
            wave += gaussian(arr[iter], gamma)
        wave = wave/tlen
        mind = np.argmax(wave)
        if mind < 100:
            snl += 1
        else:
            snr += 1

    alll = sum(his[0]) + snl
    allr = sum(his[1]) + snr
    his[0].pop(0)
    his[0].append(snl)
    his[1].pop(0)
    his[1].append(snr)
    if alll+allr > 0:
        srate[liter-8] = alll/(alll+allr)


data = np.zeros((2,slim-llim-7), dtype=np.float)
data[0] = np.linspace(7, slim-llim, slim-llim-7)
data[1] = srate

plt.figure(figsize=(10, 8))
plt.scatter(np.linspace(1, slim-llim+1, slim-llim-7), srate, s=10, marker='o')
plt.show()