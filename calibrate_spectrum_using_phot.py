#! //anaconda/bin/python

"""*******************************************************

******************************************************"""
#print __doc__


import os
import numpy as np
import pdb
import pylab
import logging
import pyphot
import get_filter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def calibrate_spectrum_using_single_phot(Spectrum,Photometry, comments='#',WaveUnits=None,lib=None,show_plot=False,name=None,filters_directory=None):
	"""Description:.
	Input  :- Either a numpy array or a path to file with 1st column: wavelemgth, anf second column flux in [WaveUnits,erg/cm^2/s/A]
            - Photometry [filter family, filter name, value of the flux]
            - WaveUnits units of wavelength in Spectrum
				- 'A' angstrom
				- 'nm' nanometers
				- 'cm' centimeters
				- 'm' meters
				- 'Hz' TO DO
				- 'eV' TO DO
				- 'keV' TO DO
			- in case Spectrum is a file, the lines starting with the specified comment symbol will be skipped
			- in case the filter through which you have the photometry is not in the PyPhot library or is not P48 r. Default is the pyphot library
			- show_plot: if set to True, show the plot. create it in anycase.
			- name: optional, name to be added to the output file, if you run this function in a loop and do not want the output directories to be removed
	Output :- 2N numpy array, with wavelength (same units as in Spectrum), and calibrated Spectrum (alpha*Spectrum)
	Tested : ?
	    By : Maayane T. Soumagnac Nov 2016
	   URL :
	Example: [Xi_array,best_temp,index_min,best_coeff,best_radius]=fit_black_body.fit_black_body_spec(black_body_with_errors,distance=1.,show_plot=True)
	Reliable:
	 TO DO: give an option for the speed: decrease the length of TempVec, give the spec units as options,
	 and change the location /Users/maayanesoumagnac/Maayane_Astro_python_library/data to the user one"""
	if isinstance(Spectrum, str):
		print('you provided a filename for the spectrum')
		Spectrum = np.genfromtxt(Spectrum, comments=comments)
	elif isinstance(Spectrum, np.ndarray):
		print('you provided a numpy array for the spectrum')
		Spectrum = Spectrum
	else:
		print('ERROR: Unsupported type for Spec')
		pdb.set_trace()
	if isinstance(Spectrum, str):
		name = '_'+os.path.splitext(os.path.basename(Spectrum))[0]
	else:
		name = ''
	#if name==None:
	if os.path.exists('./outputs_from_calibrate_spectrum_with_phot'+name):
			logger.info('output_path/txt_files did exist')#, I am removing it and creating a new one')
			#shutil.rmtree('./outputs_from_calibrate_spectrum_with_phot')
	else:
			logger.info('the output file file did not exist yet. I am creating it now')
			os.makedirs('./outputs_from_calibrate_spectrum_with_phot'+name)
	output_file='./outputs_from_calibrate_spectrum_with_phot'+name
	#else:
	#	if os.path.exists('./outputs_from_fit_black_body_flux_filters_function/'+str(name)):
	#		logger.info('output_path/txt_files did exist, I am removing it and creating a new one')
	#		shutil.rmtree('./outputs_from_fit_black_body_flux_filters_function/'+str(name))
	#	else:
	##		logger.info('the output file file did not exist yet. I am creating it now')
	#	os.makedirs('./outputs_from_fit_black_body_flux_filters_function/'+str(name))
	#	output_file='./outputs_from_fit_black_body_flux_filters_function/'+str(name)

	#******************** converts spectrum waveunits to AA ***********************
	Spectrum_in_meters=np.zeros(np.shape(Spectrum))
	if WaveUnits == None:
		print('you need to provide a WaveUnits variable')
		pdb.set_trace()
	if WaveUnits.lower()=='a':
		Spectrum_in_meters[:,0]=Spectrum[:,0]*1e-10 #A in meters
	elif WaveUnits.lower()=='nm':
		Spectrum_in_meters[:,0]=Spectrum[:,0]*1e-9 #nm in meters
	elif WaveUnits.lower()=='cm':
		Spectrum_in_meters[:,0]=Spectrum[:,0]*1e-2 #cm in meters
	elif WaveUnits.lower()=='m':
		Spectrum_in_meters[:,0]=Spectrum[:,0]
	elif WaveUnits==None:
		Spectrum_in_meters[:,0]=Spectrum[:,0] # default is meters
	else:
		print('WaveUnits unsupported yet')
		pdb.set_trace()
	wavelengths_AA=Spectrum_in_meters[:,0]*1e10

	#******************** Apply filter to spectrum ***********************

	if lib == None:
		lib = pyphot.get_library()
		print("Library contains: ", len(lib), " filters")
	else:
		print('the library was defined outside')
		lib = lib
	filter_family=Photometry[0]
	filter_name=Photometry[1]
	# print j
	#if filter_family.lower() == 'ptf_p48':
	#	if filter_name.lower() == 'r':
	#		print('You gave the R filter of the PTF_P48 family')
	#		Transmission = np.genfromtxt(
	#			'/Users/maayanesoumagnac/Maayane_Astro_python_library/data/filters/P48_R_T.rtf', delimiter=None)
	#		P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_R', dtype='photon', unit='Angstrom')
	#else:
	#	f = lib.find(filter_family.lower())  # filter family
	#		# for name in f:
	#		#    lib[name].info(show_zeropoints=True)
	#	P = lib[filter_name]  # filter name
	#Filter_vector = np.empty([1, 2], dtype=object)
	Filter_vector = [[filter_family, filter_name]]
	#Filter_vector[1] = [str('ptf_p48'), str('r')]
	[P_dic, wav] = get_filter.make_filter_object(Filter_vector,filters_directory=filters_directory)
	P=P_dic['filter_object'][0]
	print(P)
	# this calculates \bar{f}(P) through through filter P. wavelength must be given in the same units as wavelength units defined in filters
	### print('wavelengths_AA',wavelengths_AA)
	### print('Spectrum[:, 1]',Spectrum[:, 1])
	flux_through_filter = P.get_flux(wavelengths_AA, Spectrum[:, 1])
	central_wavelength=P.cl.item()
	#print('the photometry point I am using to calibrate is',Photometry[2])
	#pdb.set_trace()
	alpha = Photometry[2]/flux_through_filter
	print('The filter needs to be scaled up by a factor of {0} in order to math the photometry'.format(alpha))
	Calibrated_spectrum=np.zeros(np.shape(Spectrum))
	Calibrated_spectrum[:,0]=Spectrum[:,0]
	Calibrated_spectrum[:,1]=alpha*Spectrum[:,1]

	#CHECK
	calibrated_spectrum_flux_through_filter=P.get_flux(wavelengths_AA, Calibrated_spectrum[:, 1])
	#print calibrated_spectrum_flux_through_filter
	#print  Photometry[2]
	#print (calibrated_spectrum_flux_through_filter==Photometry[2])
	logger.info('calibrated_spectrum_flux_through_filter-Photometry[2]={0}'.format(calibrated_spectrum_flux_through_filter-Photometry[2]))
	if calibrated_spectrum_flux_through_filter-Photometry[2] >=0.00001e-16:
		logger.info('ERROR:the calibrated flux synthetic photometry does not match the photometry used to calibrate!')
		pdb.set_trace()
	#print Photometry[2]
	#pdb.set_trace()
	pylab.figure()
	pylab.plot(wavelengths_AA,Spectrum[:,1],'b',label='Uncalibrated spectrum')
	pylab.plot(wavelengths_AA,Calibrated_spectrum[:,1],'r',label='Calibrated spectrum')
	pylab.plot(central_wavelength,Photometry[2],'g*',label='Photometry used for calibration',markersize=10)
	pylab.grid()
	pylab.legend()
	pylab.savefig('./outputs_from_calibrate_spectrum_with_phot/calibration_results.pdf', facecolor='w',
				  edgecolor='w',
				  orientation='portrait', papertype=None, format='pdf', transparent=False, bbox_inches=None,
				  pad_inches=0.1)

	if show_plot==True:
		pylab.show()
	return Calibrated_spectrum