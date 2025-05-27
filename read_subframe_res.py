import numpy as n
# ncdump -h
import matplotlib.pyplot as plt

import xarray as xr
# Load the NetCDF file
# bullet at rifle range
#ds = xr.open_dataset("2024-03-11_12.nc",auto_complex=True)
# perpendicular radar failed
#ds = xr.open_dataset("2025-02-24T15:52:51_6606_perpendicular.nc",auto_complex=True)
ds = xr.open_dataset("/data2/space_radar/2025-02-25T15:39:32_6607_perpendicular.nc",auto_complex=True)
 
 
#
# number of range-gates
# n = 128
# chirp = 4096
# frame = 142 (each frame has chirp chirps)
#
# radar_cube
# maximum intermediate frequency =
# sampling frequency = sample-rate
rc=ds["radar_cube"][:]
nchirp=rc.shape[2]
nrg=rc.shape[4]

fi=ds["frame_of_interest"]

print(rc.shape)
print(ds)
print(rc.values)
subframe=256
step=32
idx=n.arange(subframe,dtype=int)
nsub=int((nchirp-subframe)/step)
imgi=0
for i in range(33,rc.shape[1]):
    deco=n.zeros([nchirp,nrg],dtype=n.complex64)
    for j in range(nchirp):
        deco[j,:]=n.fft.fftshift(n.fft.fft(rc.values[0,i,j,0,:]))
    
    for sfi in range(nsub):
        decord=n.zeros([subframe,nrg],dtype=n.complex64)
        for j in range(nrg):
            decodc=deco[idx+sfi*step,j]-n.mean(deco[idx+sfi*step,j])
            print(decodc.shape)
            decord[:,j]=n.fft.fftshift(n.fft.fft(decodc))

        nfloor=n.median(10.0*n.log10(n.abs(decord.T)**2.0))
        plt.title("frame=%d,%d"%(i,sfi))
        plt.pcolormesh(10.0*n.log10(n.abs(decord.T)**2.0),vmin=nfloor)
        plt.xlabel("Doppler (samples)")
        plt.ylabel("Range (samples)")
        plt.colorbar()
#        plt.show()
        plt.savefig("frame-%06d.png"%(imgi))
        plt.close()
        imgi+=1
#    plt.show()
        
