import numpy as np
from scipy.fftpack import fft, fftfreq, fft2
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.colors import LogNorm, Normalize, PowerNorm

nx = 128
Dx = 0.1

Dt = 0.002

lx=nx*Dx
m = 4
kt=2*np.pi/lx * m
w = (kt**2 /2) *(np.sqrt(1+4/kt**2) + 1)
finalTime = 2*np.pi / w

column_names = ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'P']

def read_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=column_names)
    return df

results_dir = 'whislerwaveres'

quantities = {
    'rho': [],
    'vx': [],
    'vy': [],
    'vz': [],
    'Bx': [],
    'By': [],
    'Bz': [],
    'P': []
}

for filename in os.listdir(results_dir):
    if filename.startswith("PRK2_") and filename.endswith(".txt"):
        df = read_file(os.path.join(results_dir, filename))

        for quantity in quantities.keys():
            quantities[quantity].append(df[quantity].values.reshape(nx))

for quantity in quantities.keys():
    quantities[quantity] = np.array(quantities[quantity])

def whistler(k):
    return (k**2 /2) * (np.sqrt(1+4/k**2) + 1)

def leftAlfvenIonCyclotron(k):
    return (k**2 /2) * (np.sqrt(1+4/k**2) - 1)

k = 2*np.pi*fftfreq(nx, Dx)
w = 2*np.pi*fftfreq(int(finalTime/Dt) + 1, Dt)

halfk = int(len(k)/2)
halfw = int(len(w)/2)
kplus = k[:halfk]
wplus = w[:halfw]

modes = [m]
kmodes = 2*np.pi*np.asarray([1/(nx*Dx) * m for m in modes])

B_combined = quantities['By'] + 1j * quantities['Bz']

# window_space = np.hanning(nx)
# window_time = np.hanning(int(finalTime / Dt) + 1)
# window = np.outer(window_time, window_space)

# B_combined *= window

B_fft = fft2(B_combined)
B_fft_abs = np.abs(B_fft)[:halfw,:halfk]

Xx, Yy = np.meshgrid(kplus, wplus)

fig, ax = plt.subplots()
pcm = ax.pcolormesh(Xx, Yy, B_fft_abs, cmap='plasma', shading='auto', norm=PowerNorm(gamma=0.3, vmin=B_fft_abs.min(), vmax=B_fft_abs.max()))
fig.colorbar(pcm, ax=ax, label='Magnitude of FFT(By + i*Bz)')
ax.plot(kplus, whistler(kplus), marker = '+', color='k')
ax.plot(kplus, leftAlfvenIonCyclotron(kplus), marker = '*', color='g')
ax.plot(kplus, kplus)
for km in kmodes:
    ax.axvline(km, ls='--', color='k')
    ax.plot(km, whistler(km), marker='o')

plt.show()