import os, sys
sys.path.append(os.path.join('G:\\Ken_RB\\GCS_Scripts', "pyhht-dev"))
from pyhht.emd import EMD
from DEM_Detrending import *

coords = station_coords('G:\\Ken_RB\\c01\\centerline\\clipped_centerline.shp', 'G:\\Ken_RB\\c01\\stationing\\stations_100.shp', 'G:\\Ken_RB\\c01\\DEM\\detrending\\EMD\\RB_DEM.tif')

df = pd.read_excel(coords)
loc = np.array(df['dist_down'])
profile = np.array(df['RASTERVALU'])
decomposer = EMD(profile)
imfs = decomposer.decompose()

print(np.shape(imfs))

fig, ax = plt.subplots(4, 1, sharex=True)

ax[0].set(title='Profile')
ax[0].plot(loc, profile, label='profile')
ax[0].plot(loc, profile-sum(imfs[:-1]), label='trend = profile - sum(imfs[:-1])')
ax[0].legend()

ax[1].set(title='IMF 1')
ax[1].plot(loc, imfs[0])

ax[2].set(title='IMF 2')
ax[2].plot(loc, imfs[1])

ax[3].set(title='detrended profile = sum(imfs[:-1])')
ax[3].plot(loc, sum(imfs[:-1]))

for i in range(4):
    ax[i].grid()
    if i != 0:
        ax[i].axhline(y=0, linestyle='--', color='black')

# plt.show()

# trend fitting
coords = df[['dist_down', 'POINT_X', 'POINT_Y', 'RASTERVALU']]
z_fits = imfs[-1][::-1]
# rename columns
coords.columns = ['distance downstream', 'x', 'y', 'z']
# get spacing between stations
spacing = abs(coords['distance downstream'][1] - coords['distance downstream'][0])
# index and sort values going downstream
coords.index = map(int, coords['distance downstream'] / spacing)
coords = coords.sort_values(by=['distance downstream'])

# station is numbered going downstream, regardless if ET GeoWizards numbered stations going upstream
station = coords.index.tolist()
# distance downstream
dist = coords['distance downstream'].tolist()
# coordinates
x = coords.x.tolist()
y = coords.y.tolist()
z = coords.z.tolist()
z_res = [z_val - z_fit for z_val, z_fit in zip(z, z_fits)]

xyz_fit = pd.DataFrame({'x': x, 'y': y, 'z_fit': z_fits})
xyz_fit.to_csv('G:\\Ken_RB\\c01\\DEM\\detrending\\EMD\\xyz_fit.csv', index=False)
detrended_DEM = detrend_DEM('G:\\Ken_RB\\c01\\DEM\\detrending\\EMD\\xyz_fit.csv', 'G:\\Ken_RB\\c01\\DEM\\detrending\\EMD\\RB_DEM.tif')
# test bivariate EMD?

matrix = arcpy.RasterToNumPyArray('G:/Ken_RB/Scratch/Rainbow_Basin_DEM2.tif')

