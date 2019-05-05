"""*******************************************************
This code plots spectra for ZTF19aanpcep - a possible Type IIn
******************************************************"""
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
color_dict['ARC']='y'
color_dict['Keck1']='g'
color_dict['WHT']='M'
z = 0.0314      # checked
linesx = {r'$H_{\alpha}$': 6563, r'$H_{\beta}$': 4861, r'$H_{\gamma}$': 4341}
    #, r'$H_{\delta}$': 4102,r'$H_{\epsilon}$': 3970, r'$H_{\zeta}$': 3889, r'$H_{\eta}$': 3835}
distance_modulus = 35.64        #checked
distance_pc = distances_conversions.DM_to_pc(distance_modulus)   # sn distance in pc
EBV = 0.014631      # checked
explosion_date = 2458568.51696456 # checked
data_dicts = Read_data.read_data_Marshall_simple(
    '/home/jinlng/test_dir/Type_IIn/ZTF19aanpcep/data_Marshal.txt',
    filters_directory=filters_directory)

P48_r_array = np.zeros((np.shape(data_dicts['jd'][data_dicts['filter']=='r_p48'])[0],2))
P48_r_array[:,0] = data_dicts['jd'][data_dicts['filter']=='r_p48']
P48_r_array[:,1] = data_dicts['flux'][data_dicts['filter']=='r_p48']


### Calibration ###
### 1 ###
spec_1 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aanpcep/ZTF19aanpcep_20190329_P60_v1.ascii',   # ok
    instrument='P60',skiprows=0,time='2019-03-29T21:53:51.729301',time_format='utc',show=False) #ok
spec_1_cal = spec_1.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)

### 2 ###
spec_2 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aanpcep/ZTF19aanpcep_20190330_APO_v1.ascii',
    instrument='ARC',skiprows=0,time='2019-03-30T07:40:52.608',time_format='utc',show=False)    #ok
spec_2_cal = spec_2.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)

### 3 ###
spec_3 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aanpcep/ZTF19aanpcep_20190330_P60_v1.ascii',   #ok
    instrument='P60',skiprows=0,time='2019-03-30T22:12:22.096502',time_format='utc',show=False)    #ok
spec_3_cal = spec_3.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)

### 4 ###
spec_4 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aanpcep/ZTF19aanpcep_20190330_P60_v2.ascii',   #ok
    instrument='P60',skiprows=0,time='2019-03-30T23:15:09.712795',time_format='utc',show=False)    #ok
spec_4_cal = spec_4.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)

### 5 ###
spec_5 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aanpcep/ZTF19aanpcep_20190331_P60_v1.ascii',   #ok
    instrument='P60',skiprows=0,time='2019-03-31T02:25:56.486801',time_format='utc',show=False)    #ok
spec_5_cal = spec_5.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)

### 6 ###
spec_6 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aanpcep/ZTF19aanpcep_20190403_Keck1_v1.ascii',
    instrument='Keck1',skiprows=0,time='2019-04-03T12:18:26',time_format='utc',show=False)    #ok
spec_6_cal = spec_6.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)


### 7 ###
spec_7 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aanpcep/ZTF19aanpcep_20190407_APO_v1.ascii',
    instrument='ARC',skiprows=0,time='2019-04-07T09:21:27.728',time_format='utc',show=False)    #ok
spec_7_cal = spec_7.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)

### 8 ###


spectra = [spec_1, spec_2, spec_3, spec_4, spec_5, spec_6, spec_7]
spectra_cal = [spec_1_cal, spec_2_cal, spec_3_cal, spec_4_cal, spec_5_cal, spec_6_cal, spec_7_cal]



dates = []
for i in spectra[::-1]:
    print('the date is',i.date_utc_hms())
    dates.append(i.date_utc_hms())

print('there are {0} spectra_cal in total'.format(len(spectra_cal)))

#offsets = np.zeros(len(spectra))
#for i in range(len(spectra_cal)):
#    offsets[i] = np.sqrt(i)*np.max(spectra_cal[0][:,1])

offsets = [0, 0.45*np.max(spectra_cal[-1][:,1]),0.45*np.max(spectra_cal[-2][:,1]),
    1.5*np.max(spectra_cal[-3][:,1]),2*np.max(spectra_cal[-4][:,1]),
    3*np.max(spectra_cal[-5][:,1]),0.25*np.max(spectra_cal[-6][:,1])]

phases = ['+'+str(round(spectra[i].date_jd()-explosion_date,2)) for i in range(len(spectra))]
annotations = [(spectra_cal[0][-1, 0]+100,np.min(spectra_cal[0][:,1])+offsets[-1]),
    (spectra_cal[1][-1, 0]+100,np.mean(spectra_cal[1][:,1])+offsets[-2]),
    (spectra_cal[2][-1, 0]+100,np.min(spectra_cal[2][:,1])+offsets[-3]),
    (spectra_cal[3][-1, 0]+100,np.min(spectra_cal[3][:,1])+offsets[-4]),
    (spectra_cal[4][-1, 0]+100,np.min(spectra_cal[4][:,1])+offsets[-5]),
    (spectra_cal[5][-1, 0]-300,np.min(spectra_cal[5][:,1])+offsets[-6]+0.1e-15),
    (spectra_cal[6][-1, 0]-100,np.min(spectra_cal[6][:,1])+offsets[-7])]
### smoothing ###
signal_smooth = []
for i in range(len(spectra_cal)):
    signal_smooth.append(signal.savgol_filter(spectra_cal[i][:,1],53,7))


#plot all spectrum
fig = plt.figure(figsize=(6.5,10.0))
for i,speci in enumerate(spectra_cal):
    print('i is',i)
    print('offsets[i] is',offsets[i])
    print('spectra[i].date_utc_hms()',spectra[i].date_utc_hms())
    print('spectra[i].instrument',spectra[i].instrument)
    plt.plot(speci[:, 0], speci[:, 1] + offsets[::-1][i],
               label=spectra[i].instrument + ', ' + spectra[i].date_utc_hms(), color=color_dict[spectra[i].instrument], alpha=0.7)
# plot each smoothed spectrum
for i,speci in enumerate(signal_smooth):
    plt.plot(spectra_cal[i][:, 0], speci + offsets[::-1][i], '-k', alpha=0.7)

for i,j in linesx.items():
    plt.axvline(j*(z+1),linestyle='-.',color='grey')
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
for i,j in enumerate(annotations):
    ax.annotate(phases[i],xy = (annotations[i]), fontsize=12.5)
plt.title('ZTF19aanpcep', fontsize=20)
plt.ylim(-0.15e-15,2.75e-15)
plt.xlim(3000,11000)
plt.ylabel(r'$\rm{flux [erg/sec/cm^2/\AA ]}$', fontsize=20)
plt.xlabel(r'wavelength ($\AA$)', fontsize=20)
plt.grid()
plt.tight_layout()
pylab.savefig('spectra_ZTF19aanpcep.pdf', figsize=(8,12), facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='pdf',transparent=False, bbox_inches='tight')
plt.show()