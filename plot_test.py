"""*******************************************************
This code plots spectra for ZTF19aaadwfi
it serves as an example for future use
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

filters_directory='/home/jinlng/test_dir/Type_IIn/Filters/Filters'
z = 0.222
linesx = {r'$H_{\alpha}$': 6563, r'$H_{\beta}$': 4861, r'$H_{\gamma}$': 4341}
    #, r'$H_{\delta}$': 4102,r'$H_{\epsilon}$': 3970, r'$H_{\zeta}$': 3889, r'$H_{\eta}$': 3835}
distance_modulus = 40.23
distance_pc = distances_conversions.DM_to_pc(distance_modulus)   # sn distance in pc
EBV = 0.012041
explosion_date = 2458478.505622309 # time of explosion given
data_dicts = Read_data.read_data_Marshall_simple(
    '/home/jinlng/test_dir/Type_IIn/ZTF19aaadwfi/data_Marshal_latest.txt',
    filters_directory=filters_directory)

P48_r_array = np.zeros((np.shape(data_dicts['jd'][data_dicts['filter']=='r_p48'])[0],2))
P48_r_array[:,0] = data_dicts['jd'][data_dicts['filter']=='r_p48']
P48_r_array[:,1] = data_dicts['flux'][data_dicts['filter']=='r_p48']


### Calibration ###
### 1 ###
Spectrum_29012019_LT = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aaadwfi/ZTF19aaadwfi_20190129_LT_v2.ascii',
    instrument='LT',skiprows=0,time='2019-01-29T03:48:37.924',time_format='utc',show=False)
Spectrum_29012019_LT_cal = Spectrum_29012019_LT.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
### 2 ###
Spectrum_04022019_T1m = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aaadwfi/ZTF19aaadwfi_20190204_LCOGT1m_v2.ascii',
    instrument='LCOGT1M',skiprows=0,time='2019-02-04T12:10:23',time_format='utc',show=False)
Spectrum_04022019_T1m_cal = Spectrum_04022019_T1m.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
### 3 ###
Spectrum_04022019_LT = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aaadwfi/ZTF19aaadwfi_20190204_LT_v1.ascii',
    instrument='LT',skiprows=0,time='2019-02-04T02:19:18.367',time_format='utc',show=False)
Spectrum_04022019_LT_cal = Spectrum_04022019_LT.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
### 4 ###
Spectrum_12022019_P200 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aaadwfi/ZTF19aaadwfi_20190212_P200_v1.ascii',
    instrument='P200',skiprows=0,time='2019-02-12T06:49:17.094',time_format='utc',show=False) 
Spectrum_12022019_P200_cal = Spectrum_12022019_P200.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
### 5 ### problematique at the high wavelength end, mask with ylim
Spectrum_12022019_P200_bis = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aaadwfi/ZTF19aaadwfi_20190212_P200_v2.ascii',
    instrument='P200',skiprows=0,time='2019-02-12T17:04:30',time_format='utc',show=False) #time not clear
Spectrum_12022019_P200_bis_cal = Spectrum_12022019_P200_bis.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
### 6 ###
Spectrum_20190315_NOT = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aaadwfi/ZTF19aaadwfi_20190315_NOT_v1.ascii',
    instrument='NOT',skiprows=0,time='2019-03-15T04:12:04.638',time_format='utc',show=False)
Spectrum_20190315_NOT_cal = Spectrum_20190315_NOT.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
### 7 ###
Spectrum_20190315_P60 = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aaadwfi/ZTF19aaadwfi_20190315_P60_v1.ascii',
    instrument='P60',skiprows=0,time='2019-03-15T01:04:59.461551',time_format='utc',show=False)
Spectrum_20190315_P60_cal = Spectrum_20190315_P60.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)



spectra = [Spectrum_29012019_LT, Spectrum_04022019_T1m, Spectrum_04022019_LT, Spectrum_12022019_P200,
            Spectrum_12022019_P200_bis, Spectrum_20190315_NOT, Spectrum_20190315_P60]
spectra_cal = [Spectrum_29012019_LT_cal, Spectrum_04022019_T1m_cal, Spectrum_04022019_LT_cal, Spectrum_12022019_P200_cal,
            Spectrum_12022019_P200_bis_cal, Spectrum_20190315_NOT_cal, Spectrum_20190315_P60_cal]

dates = []
for i in spectra[::-1]:
    print('the date is',i.date_utc_hms())
    dates.append(i.date_utc_hms())

print('there are {0} spectra_cal in total'.format(len(spectra_cal)))

offsets = np.zeros(len(spectra))
for i in range(len(spectra_cal)):
    offsets[i] = i

phases = ['+'+str(round(spectra[i].date_jd()-explosion_date,2)) for i in range(len(spectra))]
annotations = [(8100,-10.5),(10100,-12),(8100,-12.5),(10100,-13.8),(10100,-14.5),(9500,-15.2),(9300,-16.3)]

### smoothing ###
signal_smooth = []
for i in range(len(spectra_cal)):
    signal_smooth.append(signal.savgol_filter(spectra_cal[i][:,1],53,7))


fig = pylab.figure()
for i,speci in enumerate(spectra_cal):
    print('i is',i)
    print('offsets[i] is',offsets[i])
    print('spectra[i].date_utc_hms()',spectra[i].date_utc_hms())
    print('spectra[i].instrument',spectra[i].instrument)
    pylab.plot(speci[:, 0], np.log10(speci[:, 1]) + offsets[::-1][i],
               label=spectra[i].instrument + ', ' + spectra[i].date_utc_hms(), alpha=0.7)


for i,speci in enumerate(signal_smooth):
    pylab.plot(spectra_cal[i][:, 0], np.log10(speci) + offsets[::-1][i], '-k', alpha=0.7)


for i,j in linesx.items():
    pylab.axvline(j*(z+1),linestyle='-.',color='grey')
ax = pylab.gca()
handles, labels = ax.get_legend_handles_labels()
for i,j in enumerate(annotations):
    ax.annotate(phases[i],xy = (annotations[i]))
pylab.title('ZTF19aaadwfi')
pylab.ylim(-16.5,-9.7)
pylab.xlim(3000,11000)
pylab.ylabel(r'$\rm{log(flux [erg/sec/cm^2/\AA ])}$')
pylab.xlabel(r'wavelength ($\AA$)')
pylab.grid()
pylab.show()




### IF ONLY ONE SPECTRA ###
'''
### Calibration ###
### the only spectra ###
spectra = class_spectrum.spectrum(path_to_data=
    '/home/jinlng/test_dir/Type_IIn/ZTF19aadgimr_20190204_NOT_v1.ascii',    # good to go but time unclear
    instrument='NOT',skiprows=0,time='2019-02-04T00:00:00',time_format='utc',show=False)
spectra_cal = spectra.calibrate(array_of_calibrating_data=P48_r_array,
    filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)


print('there are {0} spectra_cal in total'.format(len(spectra_cal)))


phases = '+'+str(round(spectra.date_jd()-explosion_date,2)) 
annotations = [(9800,-15.5)]

### smoothing ###
signal_smooth = signal.savgol_filter(spectra_cal[:,1],53,7)

#plot spectra
fig = plt.figure()
plt.plot(spectra_cal[:, 0], np.log10(spectra_cal[:, 1]),
               label=spectra.instrument + ', ' + spectra.date_utc_hms(), alpha=0.7)
# plot smoothed spectra
plt.plot(spectra_cal[:, 0], np.log10(signal_smooth), '-k', alpha=0.7)

for i,j in linesx.items():
    plt.axvline(j*(z+1),linestyle='-.',color='grey')
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
for i,j in enumerate(annotations):
    ax.annotate(phases,xy = (annotations[i]))
plt.title('ZTF19aadgimr')
plt.ylim(-15.8,-14.6)
plt.xlim(3000,11000)
plt.ylabel(r'$\rm{log(flux [erg/sec/cm^2/\AA ])}$')
plt.xlabel(r'wavelength ($\AA$)')
plt.grid()
plt.show()

'''