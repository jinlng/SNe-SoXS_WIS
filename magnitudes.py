#! //anaconda/bin/python

"""*******************************************************
This module calcultes all kind of magnitudes
*****************************************************
"""
#print __doc__

import pdb
import numpy as np
import get_filter
import logging
import pyphot


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


##################### apparent mag ######################

def apparent_mag_to_absolute_mag(apparent_magnitude,distance=None,distance_modulus=None):
    """Description: given an array of apparent magnitude in a given filter, and a distance (pc) (or a distance modulus), calculates the array of abslute
        Input  :- an array of apparent magnitudes
                - distance in pc
        Output : array of absolute magnitude
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example:
        Reliable:  """
    absolute_magnitude=np.zeros(np.shape(apparent_magnitude))
    if distance is not None:
        #absolute_magnitude=-5*np.log10(distance)+5.+apparent_magnitude
        absolute_magnitude = apparent_magnitude-5 * np.log10(distance) + 5
    elif distance_modulus is not None:
        absolute_magnitude=apparent_magnitude-distance_modulus
    else:
        logger.info('ERROR: you need to give either the distance modulus or the distance in pc for the function to work')
        pdb.set_trace()
    return absolute_magnitude

def magAB_in_filter_to_flux_in_filter(mag_in_filter,Filter_vector=None,Filter_object=None,filters_directory=None):
    """Description: given an array of apparent magnitude_AB in a given filter, calculates the array of fluxes, ith the equation
        m_AB(P)=-2.5logf(P)-2.5log(lambda_P^2/c)-48.6
        Input  :- an array of apparent magnitudes
                - either
                    Filter vector: N-long array of of arrays [[filter family, filtername],[filter family, filtername],[filter family, filtername],...etc]
            the family and names are the ones provided by the pyphot library
                    or Filter object: a filter defined as in pyphot, e.g. P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_R', dtype='photon',
                  unit='Angstrom')
        Output : array of fluxes in erg(s AA cm^2)
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example:magnitudes.magAB_in_filter_to_flux_in_filter(mag_value,Filter_vector=np.array([['sdss','SDSS_g']]))
        Reliable:  """
    fluxes=np.zeros(np.shape(mag_in_filter))
    tens = 10.0 * np.ones(np.shape(mag_in_filter))
    if Filter_vector is not None:
        for i,j in enumerate(Filter_vector):
            print(j)
            [P,wav]=get_filter.make_filter_object([j],filters_directory=filters_directory)
            print(P['filter_object'])
            print("P['filter_object'][0].AB_zero_flux is",P['filter_object'][0].AB_zero_flux)
            fluxes=P['filter_object'][0].AB_zero_flux*np.power(tens,-mag_in_filter/2.5)
    elif Filter_object is not None:
        if isinstance(Filter_object, (list,)) and len(Filter_object)==1:
            print('Filter_object[0].AB_zero_flux is',Filter_object[0].AB_zero_flux)
            fluxes = Filter_object[0].AB_zero_flux * np.power(tens, -mag_in_filter / 2.5)
        else:
            print('Filter_object.AB_zero_flux is', Filter_object.AB_zero_flux)
            fluxes = Filter_object.AB_zero_flux * np.power(tens, -mag_in_filter / 2.5)
    else:
        print('ERROR: you need to define a filter, either a Filter_vector or a filter object as in Pyphot')
        pdb.set_trace()
    #print(fluxes)
    #print(type(fluxes))
    #print(fluxes.magnitude)
    #pdb.set_trace()
    return fluxes

def error_on_mag_into_error_on_flux(magerr,flux):
    if isinstance(magerr,np.ndarray):
        fluxerr=np.abs(-2.303/2.5*magerr*flux) ##1.0855405992184108=2.5/2.303 car ln=2.303log
    else:
        fluxerr = abs(-2.303 / 2.5 * magerr * flux)
    return fluxerr
##################### absolute mag ######################

def absolute_mag_to_apparent_mag(absolute_magnitude,distance=None,distance_modulus=None):
    """Description: given an array of absolute magnitude in a given filter, and a distance (pc) (or a distance modulus), calculates the array of apparent
        Input  :- an array of apparent magnitudes
                - distance in pc
        Output : array of apparent magnitude
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example:
        Reliable:  """
    apparent_magnitude=np.zeros(np.shape(absolute_magnitude))
    if distance is not None:
        apparent_magnitude=absolute_magnitude+5*np.log10(distance)-5
    elif distance_modulus is not None:
        apparent_magnitude=absolute_magnitude+distance_modulus
    else:
        logger.info('ERROR: you need to give either the distance modulus or the distance in pc for the function to work')
        pdb.set_trace()
    return apparent_magnitude

def absolute_mag_to_apparent_flux(absolute_magnitude,distance_pc=None,Filter_vector=None,Filter_object=None):
    """Description: given an array of absolute magnitude in a given filter, and a distance (pc) (or a distance modulus), calculates the array of apparent
        Input  :- an array of absolute magnitudes
                - distance in pc
                - either
                    Filter vector: N-long array of of arrays [[filter family, filtername],[filter family, filtername],[filter family, filtername],...etc]
            the family and names are the ones provided by the pyphot library
                    or Filter object: a filter defined as in pyphot, e.g. P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_R', dtype='photon',
                  unit='Angstrom')
        Output : array of apparent magnitude
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example:
        Reliable:  """
    print(distance_pc)
    tens = 10.0 * np.ones(np.shape(absolute_magnitude))
    #fluxes=1
    apparent_mag=absolute_mag_to_apparent_mag(absolute_magnitude,distance=distance_pc)
    print(apparent_mag)
    #pdb.set_trace()
    if Filter_vector is not None:
        lib = pyphot.get_library()
        print("Library contains:", len(lib), " filters")
        for i,j in enumerate(Filter_vector):
            print(j)
            f = lib.find(j[0].lower())#filter family
            #for name in f:
            #    lib[name].info(show_zeropoints=True)
            P = lib[j[1]]#filter
            fluxes=P.AB_zero_flux*np.power(tens,-apparent_mag/2.5)
    elif Filter_object is not None:
        fluxes = Filter_object.AB_zero_flux * np.power(tens, -apparent_mag/ 2.5)
    return fluxes
