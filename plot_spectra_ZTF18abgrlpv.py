"""********************************************************
single-spectra
********************************************************"""
__author__ = 'Jingyi'

print(__doc__)

import distances_conversions
import Read_data
import class_spectrum
import numpy as np
import pylab
import pdb
import pyphot
from scipy import signal
import matplotlib.pyplot as plt

filters_directory='/home/jinlng/test_dir/Type_IIn/Filters/Filters'
color_dict={}
color_dict['NOT']='r'
color_dict['LT']='b'
color_dict['P60']='#e17701'
color_dict['P200']='#ff028d'
color_dict['LCOGT1M']='C'
color_dict['APO']='y'
color_dict['Keck1']='g'
color_dict['WHT']='M'
### ZTF18abgrlpv ###
z = 0.2104     # checked
linesx = {r'$H_{\alpha}$': 6563, r'$H_{\beta}$': 4861, r'$H_{\gamma}$': 4341}
    #, r'$H_{\delta}$': 4102,r'$H_{\epsilon}$': 3970, r'$H_{\zeta}$': 3889, r'$H_{\eta}$': 3835}
distance_modulus = 40.10       # checked
distance_pc = distances_conversions.DM_to_pc(distance_modulus)   # sn distance in pc
EBV = 0.055386      # checked
explosion_date = 2458311.479397956 # checked
data_dicts = Read_data.read_data_Marshall_simple(
    '/home/jinlng/test_dir/Type_IIn/ZTF18abgrlpv_data_Marshal.txt', # good to go ???
    filters_directory=filters_directory)

P48_r_array = np.zeros((np.shape(data_dicts['jd'][data_dicts['filter']=='r_p48'])[0],2))
P48_r_array[:,0] = data_dicts['jd'][data_dicts['filter']=='r_p48']
P48_r_array[:,1] = data_dicts['flux'][data_dicts['filter']=='r_p48']

### Calibration ###
### the only spectra ###
spectra = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF18abgrlpv_20180717_P200_v1.ascii',     # good to go
    instrument='P200',skiprows=0,time='2018-07-17T09:34:39.510',time_format='utc',show=False)    # good to go
spectra_cal = spectra.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)


print('there are {0} spectra_cal in total'.format(len(spectra_cal)))

phases = '+'+str(round(spectra.date_jd()-explosion_date,2)) 

### smoothing ###
signal_smooth = signal.savgol_filter(spectra_cal[:,1],53,7)

#plot spectra
fig = plt.figure()
plt.plot(spectra_cal[:, 0], spectra_cal[:, 1],
               label=spectra.instrument + ', ' + spectra.date_utc_hms(), color=color_dict[spectra.instrument], alpha=0.7)
# plot smoothed spectra
plt.plot(spectra_cal[:, 0], signal_smooth, '-k', alpha=0.7)

# mark the sign for telluric
plt.plot(7600,np.min(spectra_cal[:, 1])*0.5,'*b')

for i,j in linesx.items():
    plt.axvline(j*(z+1),linestyle='-.',color='grey')
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
ax.annotate(phases,xy = ((9300,0.5e-48)), fontsize=13)
plt.title('ZTF18abgrlpv', fontsize=20)
#pylab.ylim(-16.5,-9.7)
plt.xlim(3000,11000)
plt.ylabel(r'$\rm{flux [erg/sec/cm^2/\AA ]}$', fontsize=20)
plt.xlabel(r'wavelength ($\AA$)', fontsize=20)
plt.grid()
plt.tight_layout()
pylab.savefig('spectra_ZTF18abgrlpv.pdf', facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='pdf',transparent=False, bbox_inches='tight')
plt.show()
