"""********************************************************
 NOT WORKING IN TERMINAL BUT WORKS IN JUPYTER NOTEBOOK
 deleted the NaN valeur in NOT ascii, now it works fine
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
### ZTF18acwzyor ###
z = 0.02736     # checked
linesx = {r'$H_{\alpha}$': 6563, r'$H_{\beta}$': 4861, r'$H_{\gamma}$': 4341}
    #, r'$H_{\delta}$': 4102,r'$H_{\epsilon}$': 3970, r'$H_{\zeta}$': 3889, r'$H_{\eta}$': 3835}
distance_modulus = 35.33       # checked
distance_pc = distances_conversions.DM_to_pc(distance_modulus)   # sn distance in pc
EBV = 0.045197      # checked
explosion_date = 2458472.4496397697 # 
data_dicts = Read_data.read_data_Marshall_simple(
    '/home/jinlng/test_dir/Type_IIn/ZTF18acwzyor/data_Marshal.txt',
    filters_directory=filters_directory)

P48_r_array = np.zeros((np.shape(data_dicts['jd'][data_dicts['filter']=='r_p48'])[0],2))
P48_r_array[:,0] = data_dicts['jd'][data_dicts['filter']=='r_p48']
P48_r_array[:,1] = data_dicts['flux'][data_dicts['filter']=='r_p48']

### Calibration ###
### 1 ###
spec_1 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF18acwzyor/ZTF18acwzyor_20181224_P60_v1.ascii',
    instrument='P60',skiprows=0,time='2018-12-24T06:56:25.410595',time_format='utc',show=False)
spec_1_cal = spec_1.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)

### 2 ###
### time not clear
spec_2 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF18acwzyor/ZTF18acwzyor_20190119_NOT_v1.ascii',
    instrument='NOT',skiprows=0,time='2019-01-19T06:56:25.410595',time_format='utc',show=False)
spec_2_cal = spec_2.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)

### 3 ###
spec_3 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF18acwzyor/ZTF18acwzyor_20190220_WHT_v1.ascii',
    instrument='WHT',skiprows=0,time='2019-02-20T23:55:11.797',time_format='utc',show=False)
spec_3_cal = spec_3.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)

### 4 ### this is before explosion time, not usable
#spec_4 = class_spectrum.spectrum(path_to_data=
#    '/home/jinlng/test_dir/Type_IIn/ZTF18acwzyor/tns_2018kag_2018-12-18_01-09-32_LT_SPRAT.txt',   #ok
#    instrument='LT',skiprows=0,time='2018-12-18T01:09:32',time_format='utc',show=False)    #ok
#spec_4_cal = spec_4.calibrate(array_of_calibrating_data=P48_r_array,
#    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)



spectra = [spec_1, spec_2, spec_3]
spectra_cal = [spec_1_cal, spec_2_cal, spec_3_cal]

dates = []
for i in spectra[::-1]:
    print('the date is',i.date_utc_hms())
    dates.append(i.date_utc_hms())

print('there are {0} spectra_cal in total'.format(len(spectra_cal)))


offsets = [0,np.max(spectra_cal[1][:,1])*0.2,
    np.max(spectra_cal[0][:,1])*0.2+np.max(spectra_cal[1][:,1])*0.3] 
#    np.max(spectra_cal[0][:,1])*0.2+np.max(spectra_cal[1][:,1])*0.3+np.max(spectra_cal[2][:,1])*0.3]
#offsets = np.zeros(len(spectra))
#for i in range(len(spectra_cal)):
#    offsets[i] = i*np.max(spectra_cal[i-1][:,1])

phases = ['+'+str(round(spectra[i].date_jd()-explosion_date,2)) for i in range(len(spectra))]
annotations = [(spec_1_cal[-1, 0]+100,np.min(spec_1_cal)+offsets[2]),
            (spec_2_cal[-1, 0]+100,np.min(spec_2_cal)+offsets[1]),
            (spec_3_cal[-1, 0]+100,np.min(spec_3_cal)+offsets[0])]

### smoothing ###
signal_smooth = []
for i in range(len(spectra_cal)):
    signal_smooth.append(signal.savgol_filter(spectra_cal[i][:,1],53,7))


#plot all spectrum
fig = plt.figure()
for i,speci in enumerate(spectra_cal):
    print('i is',i)
    print('offsets[i] is',offsets[i])
    print('spectra[i].date_utc_hms()',spectra[i].date_utc_hms())
    print('spectra[i].instrument',spectra[i].instrument)
    plt.plot(speci[:, 0], speci[:, 1] + offsets[::-1][i],
               label=spectra[i].instrument + ', ' + spectra[i].date_utc_hms(), 
               color=color_dict[spectra[i].instrument], alpha=0.7)
# plot each smoothed spectrum
for i,speci in enumerate(signal_smooth):
    plt.plot(spectra_cal[i][:, 0], speci + offsets[::-1][i], '-k', alpha=0.7)

# mark the sign for telluric
plt.plot(7615,np.min(spec_3_cal[:,1]*0.01),'*')

for i,j in linesx.items():
    plt.axvline(j*(z+1),linestyle='-.',color='grey')
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
for i,j in enumerate(annotations):
    ax.annotate(phases[i],xy = (annotations[i]), fontsize=13)
plt.title('SN 2018kag', fontsize=20)
#plt.ylim(-17.5,-12.5)
#plt.ylim(np.min(spec_1_cal[:,1])*0.5,np.max(spec_3_cal[:,1])+offsets[2]*1.1)  
plt.xlim(3000,11000)
plt.ylabel(r'$\rm{flux [erg/sec/cm^2/\AA ]}$', fontsize=20)
plt.xlabel(r'wavelength ($\AA$)', fontsize=20)
plt.grid()
plt.tight_layout()
pylab.savefig('spectra_ZTF18acwzyor.pdf', facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='pdf',transparent=False, bbox_inches='tight')
plt.show()
