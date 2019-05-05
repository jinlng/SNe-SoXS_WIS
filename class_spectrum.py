"""*******************************************************
This code has a class spectrum
*****************************************************
"""
#print __doc__


from astropy.time import Time
from astropy.io import ascii
from scipy.interpolate import interp1d
import read_data_from_file
import calibrate_spectrum_using_phot
import pylab
import numpy as np



class spectrum(object):
    def __init__(self,instrument,path_to_data,skiprows,time,time_format,show=True):
        self.path_to_data=path_to_data
        self.instrument=instrument
        self.time=time
        self.time_format=time_format
        self.skiprows=skiprows
        self.show=show
    def date_jd(self):
        if self.time_format=='mjd':
            date_jd=self.time+2400000.5
            print('date is',date_jd)
            return date_jd
        elif self.time_format=='utc':
            t= Time(self.time, format='isot', scale='utc')
            date_jd=t.jd
            return date_jd
        elif self.time_format=='jd':
            date_jd=self.time
            print('date_jd is',date_jd)
            return date_jd
    def date_utc_hm(self):
        if self.time_format == 'mjd':
            t=Time(self.time, format='mjd', scale='utc')
            date_utc = t.iso
            date_utc_simple=Time(date_utc,in_subfmt='date_hms',out_subfmt='date_hm').iso
            #print(date_utc_simple)
            #pdb.set_trace()
            return date_utc_simple
        elif self.time_format == 'utc':
            t = Time(self.time, format='isot', scale='utc',in_subfmt='date_hms',out_subfmt='date_hm')
            #print(t.iso)
            #pdb.set_trace()
            return t.iso
        elif self.time_format == 'jd':
            t=Time(self.time, format='jd', scale='utc')
            date_utc = t.iso
            date_utc_simple=Time(date_utc,in_subfmt='date_hms',out_subfmt='date_hm').iso
            return date_utc_simple
    def date_utc_hms(self):
        if self.time_format == 'mjd':
            t=Time(self.time, format='mjd', scale='utc')
            date_utc = t.iso
            #print(date_utc_simple)
            #pdb.set_trace()
            return date_utc
        elif self.time_format == 'utc':
            t = Time(self.time, format='isot', scale='utc')
            #print(t.iso)
            #pdb.set_trace()
            return t.iso
        elif self.time_format=='jd':
            t = Time(self.time, format='jd', scale='utc')
            date_utc = t.iso
            return date_utc
    def numpy_array(self):
        a=ascii.read(self.path_to_data)
        array=np.array(list(zip(a[0][:],a[1][:])))
        #array=read_data_from_file.read_data_into_numpy_array(self.path_to_data,skiprows=self.skiprows,delimiter=' ')
        return array
    def calibrate(self,array_of_calibrating_data,filter_family,filter_name,filters_directory=None):
        jd=array_of_calibrating_data[:,0]
        flux=array_of_calibrating_data[:,1]
        interpolated_flux = interp1d(jd, flux)  #, bounds_error=False, fill_value='extrapolate') 
        photometric_flux_interpolated_over_spectrum_day=interpolated_flux(self.date_jd()) 
        calibrated_spectrum=calibrate_spectrum_using_phot.calibrate_spectrum_using_single_phot(self.numpy_array(),
            [filter_family, filter_name, photometric_flux_interpolated_over_spectrum_day], comments='#', WaveUnits='A',
            lib=None, show_plot=False, name=None,filters_directory=filters_directory)
        #pylab.figure()
        #pylab.plot(self.numpy_array()[:, 0], self.numpy_array()[:, 1],label='original')
        #pylab.plot(calibrated_spectrum[:,0],calibrated_spectrum[:,1],label='calibrated')
        #pylab.title(self.instrument + ',t(jd)=' + str(self.date_jd()) + ',t(UTC)=' + self.date_utc())
        #pylab.grid()
        #if self.show == True:
        #    pylab.show()
        return calibrated_spectrum
    def plot(self):
        pylab.figure()
        print('self.date_jd()',self.date_jd())
        print(str(self.date_jd()))
        print('self.date_utc_hms()',self.date_utc_hms())
        #pdb.set_trace()
        pylab.plot(self.numpy_array()[:,0],self.numpy_array()[:,1])
        pylab.title(self.instrument+',t(jd)='+str(self.date_jd())+',t(UTC)='+self.date_utc_hms())
        pylab.grid()
        if self.show==True:
            pylab.show()







