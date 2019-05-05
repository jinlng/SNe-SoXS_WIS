
"""*******************************************************
This code plots spectra for 2018fif
******************************************************"""
__author__ = 'maayanesoumagnac'

print(__doc__)


import distances_conversions
import Read_data
import class_spectrum
import numpy as np
import pylab

#when I will automate this, here are the arguments: Ra, Dec, distance modulus, z, EBV, datePspectra

filters_directory='/Users/maayanesoumagnac/AstroCosmo/GitHub_repositories/PhotoFit/PhotoFit/Filters'
# INFO on 2018fif
ra= 2.360644 #15h54m53.04000s   +03d32m07.5156s
dec=47.354093
z=0.017189
linesx = {r'$H_{\alpha}$': 6563, r'$H_{\beta}$': 4861, r'$H_{\gamma}$': 4341}#, r'$H_{\delta}$': 4102,r'$H_{\epsilon}$': 3970, r'$H_{\zeta}$': 3889, r'$H_{\eta}$': 3835}
#sky_level
distance_modulus=34.3
distance_pc=distances_conversions.DM_to_pc(distance_modulus)   # sn distance in pc

EBV=0.1# with Eran's sky_ebv.m from Gauss

data_dicts=Read_data.read_data_Marshall('./data/data_Marshal_latest.txt')

P48_r_array=np.zeros((np.shape(data_dicts['r']['jd'])[0],2))
P48_r_array[:,0]=data_dicts['r']['jd']
P48_r_array[:,1]=data_dicts['r']['flux']


#Calibrate all the spectra

Spectrum_21082018_Gemini=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20180821_Gemini_N_v2.ascii.txt',instrument='Gemini',skiprows=None,time=58351.5520816792,time_format='mjd',show=False)
Spectrum_21082018_Gemini_cal=Spectrum_21082018_Gemini.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
Spectrum_21082018_Gemini.plot()
#pylab.show()
Spectrum_21082018_P200=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20180821_P200_v1.ascii',instrument='P200',skiprows=2,time='2018-08-21T12:08:32.457',time_format='utc',show=False)
Spectrum_21082018_P200_cal=Spectrum_21082018_P200.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
Spectrum_21082018_P200.plot()

Spectrum_23082018_SEDM=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20180823_P60_v2.ascii',instrument='P60',skiprows=156,time=58353.2079007,time_format='mjd',show=False)
Spectrum_23082018_SEDM_cal=Spectrum_23082018_SEDM.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_23082018_SEDM.plot()

Spectrum_25082018_LT=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20180825_LT_v1.ascii',instrument='LT',skiprows=270,time='2018-08-25T23:25:40.425',time_format='utc',show=False)
Spectrum_25082018_LT_cal=Spectrum_25082018_LT.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_25082018_LT.plot()

Spectrum_27082018_SEDM=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20180827_P60_v1.ascii',instrument='P60',skiprows=158,time=58357.1821648,time_format='mjd',show=False)
Spectrum_27082018_SEDM_cal=Spectrum_27082018_SEDM.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_27082018_SEDM.plot()

Spectrum_29082018_SEDM=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20180829_P60_v1.ascii',instrument='P60',skiprows=170,time=58359.4739745,time_format='mjd',show=False)
Spectrum_29082018_SEDM_cal=Spectrum_29082018_SEDM.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_29082018_SEDM.plot()

Spectrum_04092018_NOT=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20180904_NOT_v1.ascii',instrument='NOT',skiprows=None,time='2018-09-04T00:00:00',time_format='utc',show=False)
Spectrum_04092018_NOT_cal=Spectrum_04092018_NOT.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_04092018_NOT.plot()

Spectrum_25092018_SEDM=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20180925_P60_v1.ascii',instrument='P60',skiprows=173,time=58386.35640949989,time_format='mjd',show=False)
Spectrum_25092018_SEDM_cal=Spectrum_25092018_SEDM.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_25092018_SEDM.plot()

Spectrum_03112018_SEDM=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20181103_P60_v1.ascii.txt',instrument='P60',skiprows=181,time=58425.11831830023,time_format='mjd',show=False)
Spectrum_03112018_SEDM_cal=Spectrum_03112018_SEDM.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_03112018_SEDM.plot()

Spectrum_14112018_SEDM=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20181114_P60_v1.ascii.txt',instrument='P60',skiprows=181,time=58436.32909470005,time_format='mjd',show=False)
Spectrum_14112018_SEDM_cal=Spectrum_14112018_SEDM.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_14112018_SEDM.plot()

Spectrum_19112018_SEDM=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20181119_P60_v1.ascii.txt',instrument='P60',skiprows=181,time=58441.26803589985,time_format='mjd',show=False)
Spectrum_19112018_SEDM_cal=Spectrum_19112018_SEDM.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_19112018_SEDM.plot()

Spectrum_26112018_SEDM=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20181126_P60_v1.ascii.txt',instrument='P60',skiprows=181,time=58448.19395820005,time_format='mjd',show=False)
Spectrum_26112018_SEDM_cal=Spectrum_26112018_SEDM.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_26112018_SEDM.plot()

Spectrum_04122018_SEDM=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20181204_P60_v1.ascii.txt',instrument='P60',skiprows=181,time=58456.32505029999,time_format='mjd',show=False)
Spectrum_04122018_SEDM_cal=Spectrum_04122018_SEDM.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_04122018_SEDM.plot()

Spectrum_17122018_WHT=class_spectrum.spectrum(path_to_data='./data/spectroscopy/ZTF18abokyfk/ZTF18abokyfk_20181217_WHT_v1.ascii.txt',instrument='WHT',skiprows=None,time=2458470.334540,time_format='jd',show=False)#17/12/2018/20/1/44.26  day/month/year/hour/minute/second
Spectrum_17122018_WHT_cal=Spectrum_17122018_WHT.calibrate(array_of_calibrating_data=P48_r_array,filter_family='ztf_p48',filter_name='r_p48',filters_directory=filters_directory)
#Spectrum_17122018_WHT.plot()


dates=[]
spectra=[Spectrum_21082018_P200,Spectrum_21082018_Gemini,Spectrum_23082018_SEDM,Spectrum_27082018_SEDM,Spectrum_29082018_SEDM,Spectrum_04092018_NOT,Spectrum_25092018_SEDM,Spectrum_03112018_SEDM,Spectrum_26112018_SEDM, Spectrum_04122018_SEDM,Spectrum_17122018_WHT]
#spectra_calx=[Spectrum_21082018_P200_cal,Spectrum_21082018_Gemini_cal,Spectrum_23082018_SEDM_cal,Spectrum_27082018_SEDM_cal,Spectrum_29082018_SEDM_cal,Spectrum_04092018_NOT_cal,Spectrum_25092018_SEDM_cal,Spectrum_03112018_SEDM_cal,Spectrum_14112018_SEDM_cal,Spectrum_19112018_SEDM_cal,Spectrum_26112018_SEDM_cal,Spectrum_04122018_SEDM_cal,Spectrum_17122018_WHT_cal]
spectra_cal=[Spectrum_21082018_P200_cal,Spectrum_21082018_Gemini_cal,Spectrum_23082018_SEDM_cal,Spectrum_27082018_SEDM_cal,Spectrum_29082018_SEDM_cal,Spectrum_04092018_NOT_cal,Spectrum_25092018_SEDM_cal,Spectrum_03112018_SEDM_cal,Spectrum_26112018_SEDM_cal,Spectrum_04122018_SEDM_cal,Spectrum_17122018_WHT_cal]
#print(len(spectra))
#print(len(spectra_cal))
#pdb.set_trace()

for i in spectra[::-1]:
    print('the date is',i.date_utc_hms())
    dates.append(i.date_utc_hms())
#print(dates)
#pdb.set_trace()
'''
pylab.figure()
for spec in spectra:
    pylab.plot(spec.numpy_array()[:, 0], spec.numpy_array()[:, 1],
               label=spec.instrument + ', ' + spec.date_utc_hm(), alpha=0.8)
pylab.legend()
pylab.grid()
pylab.show()
'''
#offsets=[0+s*0.15e-15 for s in range(9)]
#offsets=[0+s*0.2 for s in range(9)]
print('there are {0} spectra_cal in total'.format(len(spectra_cal)))
offsets=np.zeros(len(spectra))
offsets[0]=0
offsets[1]=0.35
offsets[2]=0.40
#offsets[3]=0.5
offsets[3]=0.6
offsets[4]=0.7
offsets[5]=0.8
offsets[6]=1.
offsets[7]=1.4
offsets[8]=1.8
offsets[9]=2.5
offsets[10]=3.
#offsets[12]=3.5

offsets=offsets+1.
print('offsets:',offsets)
#pdb.set_trace()

pylab.figure()
pylab.plot(spectra_cal[-1][:, 0],spectra_cal[-1][:,1])
pylab.plot(spectra_cal[-1][:, 0],spectra_cal[-1][:,1])
#ax=pylab.gca()
#ax.set_yscale("log", nonposy='clip')
#pylab.show()
#pdb.set_trace()

'''
fig=pylab.figure(figsize=(6,12))
for i,speci in enumerate(spectra_cal):
    print('offsets[i] is',offsets[i])
    print('spectra_cal[i].date_utc_hm()',spectra[i].date_utc_hm())
    print('spectra[i].instrument',spectra[i].instrument)
    pylab.plot(speci[:, 0], np.log10(speci[:, 1]) + offsets[i],
               label=spectra[i].instrument + ', ' + spectra[i].date_utc_hm())#, alpha=0.9)
for i,j in linesx.items():
    pylab.axvline(j*(z+1),linestyle='-.',color='grey')
ax=pylab.gca()
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper right')
pylab.ylabel(r'$\rm{log(flux [erg/sec/cm^2/\AA ])}$')
pylab.xlabel(r'wavelength ($\AA$)')
#pylab.legend(loc='upper right')
pylab.grid()
pylab.ylim(-16.5,-11)
#ax=pylab.gca()
#ax.set_yscale("log", nonposy='clip')
pylab.tight_layout()
#plt.figure()
pylab.savefig('./pdf_files/spectra_calib.pdf', facecolor='w', edgecolor='w',
              orientation='portrait',
              papertype=None, format='pdf', transparent=False, bbox_inches=None, pad_inches=0.5)
'''

explosion_date=2458351.653729907237

phases=['+'+str(round(spectra[i].date_jd()-explosion_date,2)) for i in range(len(spectra))]
print(phases)
#pdb.set_trace()
annotations=[(9250,-12.5),(7000,-12.7),(9000,-13.4),(9200,-13.8),(9300,-13.9),(9700,-14.2),(9350,-14.35),(9300,-14.6),(9300,-14.8),(9200,-15.2),(9600,-15.5)]

fig=pylab.figure(figsize=(6,8))
for i,speci in enumerate(spectra_cal):
    print('i is',i)
    print('offsets[i] is',offsets[i])
    print('spectra[i].date_utc_hm()',spectra[i].date_utc_hm())
    print('spectra[i].instrument',spectra[i].instrument)
    pylab.plot(speci[:, 0], np.log10(speci[:, 1]) + offsets[::-1][i],
               label=spectra[i].instrument + ', ' + spectra[i].date_utc_hm())#, alpha=0.9)
for i,j in linesx.items():
    pylab.axvline(j*(z+1),linestyle='-.',color='grey')
ax=pylab.gca()
handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels, loc='upper right',prop={'size': 9})
for i,j in enumerate(annotations):
    ax.annotate(phases[i],xy=(annotations[i]))
pylab.ylabel(r'$\rm{log(flux [erg/sec/cm^2/\AA ])}$')
pylab.xlabel(r'wavelength ($\AA$)')
#pylab.legend(loc='upper right')
pylab.grid()
pylab.ylim(-16,-11)
pylab.xlim(3000,11000)
#ax=pylab.gca()
#ax.set_yscale("log", nonposy='clip')
pylab.tight_layout()
#plt.figure()
pylab.savefig('./pdf_files/spectra_calib.pdf', facecolor='w', edgecolor='w',
              orientation='portrait',
              papertype=None, format='pdf', transparent=False, bbox_inches=None, pad_inches=0.5)

'''
fig=pylab.figure(figsize=(6,8))
for i,speci in enumerate(spectra_cal):
    pylab.plot(speci[:, 0], speci[:, 1],
               label=spectra[i].instrument + ', ' + spectra[i].date_utc_hm(), alpha=0.9)
#    pylab.plot(speci[:, 0], speci[:, 1] + offsets[i],
#               label=spectra[i].instrument + ', ' + spectra[i].date_utc_hm(), alpha=0.9)
for i,j in linesx.items():
    pylab.axvline(j*(z+1),linestyle='-.',color='grey')
pylab.ylabel(r'flux $[erg/sec/cm^2/\AA ]$')
pylab.xlabel(r'wavelength ($\AA$)')
pylab.legend()
pylab.grid()
ax=pylab.gca()
ax.set_yscale("log", nonposy='clip')
pylab.tight_layout()
#plt.figure()
pylab.savefig('./pdf_files/spectra_calib.pdf', facecolor='w', edgecolor='w',
              orientation='portrait',
              papertype=None, format='pdf', transparent=False, bbox_inches=None, pad_inches=0.5)
'''

pylab.show()

'''
pylab.figure()
pylab.plot(Spectrum_21082018_Gemini.numpy_array()[:,0],Spectrum_21082018_Gemini.numpy_array()[:,1],label=Spectrum_21082018_Gemini.instrument+', '+Spectrum_21082018_Gemini.date_utc_hm(),alpha=0.8)
pylab.plot(Spectrum_21082018_P200.numpy_array()[:,0],Spectrum_21082018_P200.numpy_array()[:,1],label=Spectrum_21082018_P200.instrument+', '+Spectrum_21082018_P200.date_utc_hm(),alpha=0.8)
pylab.plot(Spectrum_23082018_SEDM.numpy_array()[:,0],Spectrum_23082018_SEDM.numpy_array()[:,1],label=Spectrum_23082018_SEDM.instrument+', '+Spectrum_23082018_SEDM.date_utc_hm(),alpha=0.8)
#pylab.plot(Spectrum_25082018_LT.numpy_array()[:,0],Spectrum_25082018_LT.numpy_array()[:,1],label=Spectrum_25082018_LT.instrument,alpha=0.8)
pylab.plot(Spectrum_27082018_SEDM.numpy_array()[:,0],Spectrum_27082018_SEDM.numpy_array()[:,1],label=Spectrum_27082018_SEDM.instrument+', '+Spectrum_27082018_SEDM.date_utc_hm(),alpha=0.8)
pylab.plot(Spectrum_29082018_SEDM.numpy_array()[:,0],Spectrum_29082018_SEDM.numpy_array()[:,1],label=Spectrum_29082018_SEDM.instrument+', '+Spectrum_29082018_SEDM.date_utc_hm(),alpha=0.8)
pylab.plot(Spectrum_04092018_NOT.numpy_array()[:,0],Spectrum_04092018_NOT.numpy_array()[:,1],label=Spectrum_04092018_NOT.instrument+', '+Spectrum_04092018_NOT.date_utc_hm(),alpha=0.8)
pylab.plot(Spectrum_25092018_SEDM.numpy_array()[:,0],Spectrum_25092018_SEDM.numpy_array()[:,1],label=Spectrum_25092018_SEDM.instrument+', '+Spectrum_25092018_SEDM.date_utc_hm(),alpha=0.8)
pylab.legend()
pylab.grid()
#pylab.show()
'''

'''
offsets=[0+i*0.15e-15 for i in range(9)]
fig=pylab.figure(figsize=(6,8))
pylab.plot(Spectrum_21082018_P200_cal[:,0],Spectrum_21082018_P200_cal[:,1]+offsets[0],label=Spectrum_21082018_P200.instrument+', '+Spectrum_21082018_P200.date_utc_hm(),alpha=0.9)
pylab.plot(Spectrum_21082018_Gemini_cal[:,0],Spectrum_21082018_Gemini_cal[:,1]+offsets[1],label=Spectrum_21082018_Gemini.instrument+', '+Spectrum_21082018_Gemini.date_utc_hm(),alpha=0.9)
pylab.plot(Spectrum_23082018_SEDM_cal[:,0],Spectrum_23082018_SEDM_cal[:,1]+offsets[2],label=Spectrum_23082018_SEDM.instrument+', '+Spectrum_23082018_SEDM.date_utc_hm(),alpha=0.9)
#pylab.plot(Spectrum_25082018_LT_cal[:,0],Spectrum_25082018_LT_cal[:,1],label=Spectrum_25082018_LT.instrument,alpha=0.8)
pylab.plot(Spectrum_27082018_SEDM_cal[:,0],Spectrum_27082018_SEDM_cal[:,1]+offsets[3],label=Spectrum_27082018_SEDM.instrument+', '+Spectrum_27082018_SEDM.date_utc_hm(),alpha=0.9)
pylab.plot(Spectrum_29082018_SEDM_cal[:,0],Spectrum_29082018_SEDM_cal[:,1]+offsets[4],label=Spectrum_29082018_SEDM.instrument+', '+Spectrum_29082018_SEDM.date_utc_hm(),alpha=0.9)
pylab.plot(Spectrum_04092018_NOT_cal[:,0],Spectrum_04092018_NOT_cal[:,1]+offsets[5],label=Spectrum_04092018_NOT.instrument+', '+Spectrum_04092018_NOT.date_utc_hm(),alpha=0.9)
pylab.plot(Spectrum_25092018_SEDM_cal[:,0],Spectrum_25092018_SEDM_cal[:,1]+offsets[6],label=Spectrum_25092018_SEDM.instrument+', '+Spectrum_25092018_SEDM.date_utc_hm(),alpha=0.9)
#pylab.title('Calibrated Spectra')
for i,j in linesx.items():
    pylab.axvline(j*(z+1),linestyle='-.',color='grey')
pylab.ylabel(r'flux $[erg/sec/cm^2/\AA ]$')
pylab.xlabel(r'wavelength ($\AA$)')
pylab.legend()
pylab.grid()
pylab.tight_layout()
#plt.figure()
pylab.savefig('./pdf_files/spectra_calib.pdf', facecolor='w', edgecolor='w',
              orientation='portrait',
              papertype=None, format='pdf', transparent=False, bbox_inches=None, pad_inches=0.5)

pylab.show()
#Calibration and correction of the spectra
'''





