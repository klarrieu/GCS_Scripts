from DEM_Detrending import *
from scipy.optimize import curve_fit
import numpy as np

# coords = station_coords('G:\\Ken_RB\\c01\\centerline\\clipped_centerline.shp', 'G:\\Ken_RB\\c01\\stationing\\stations_100.shp', 'G:\\Ken_RB\\c01\\DEM\\detrending\\exp\\RB_DEM.tif')

coords = 'G:\\Ken_RB\\c01\\DEM\\detrending\\exp\\intersection_coords.xls'

df = pd.read_excel(coords)
df = df.sort_values(['dist_down'])
print(df)
loc = np.array(df['dist_down'])
profile = np.array(df['RASTERVALU'])


def func(x, a, b, c):
    return a * np.exp(- b * x) + c
    # return a * np.exp(-b * x) + c


print('fitting...')
params, cov = curve_fit(func, loc, profile, maxfev=10000)

print(params)
print(cov)

fit = [func(val, *params) for val in loc]

fig, ax = plt.subplots(2, 1, sharex=True)

ax[0].set(title='Profile')
ax[0].plot(loc, profile, label='profile')
ax[0].plot(loc, fit, label='a * exp(- b * x) + c fit')
ax[0].legend()

ax[1].plot(loc, [z_val-z_fit for z_val, z_fit in zip(profile, fit)])

for i in range(2):
    ax[i].grid()
    if i != 0:
        ax[i].axhline(y=0, linestyle='--', color='black')

plt.show()