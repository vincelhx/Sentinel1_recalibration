from lxml import objectify
import numpy as np


def get_all_geap_gains(path_aux_cal):
    """
    Find gains Geap associated with mode product and slice number from AUX_CAL.

    DOC : `https://sentinel.esa.int/documents/247904/1877131/DI-MPC-PB-0241-3-10_Sentinel-1IPFAuxiliaryProductSpecification.pdf/ae025687-c3e3-6ab0-de8d-d9cf58657431?t=1669115416469`

    Parameters
    ----------
    path_aux_cal: str 

    Returns
    ----------
    dict
        return a dict for the given (mode+pols).
        this dictionnary contains a dict with offboresight angle values and associated gains values
    """

    # Parse the XML file
    with open(path_aux_cal, 'rb') as file:
        xml = file.read()
    root_aux = objectify.fromstring(xml)

    dict_gains = {}
    for calibrationParams in root_aux.calibrationParamsList.getchildren():
        # get swath and pol
        swath = calibrationParams.swath
        polarisation = calibrationParams.polarisation

        dict_temp = {}

        increment = calibrationParams.elevationAntennaPattern.elevationAngleIncrement
        # according to the documentation, "values separated by spaces in the order I Q I Q I Q"
        valuesIQ = np.array([float(
            e) for e in calibrationParams.elevationAntennaPattern['values'].text.split(' ')])
        gain = np.sqrt(valuesIQ[::2]**2+valuesIQ[1::2]**2)
        count = gain.size

        # [-15 ; 0.05 step, 15] deg
        ang = np.linspace(-((count - 1)/2) * increment,
                          ((count - 1)/2) * increment, count)

        dict_temp["offboresightAngle"] = ang
        dict_temp["gain"] = gain
        dict_gains[swath+"_"+polarisation] = dict_temp

    return dict_gains


def get_geap_gains(path_aux_cal, mode, pols):
    """
    Find gains Geap associated with mode product and slice number from AUX_CAL.

    DOC : `https://sentinel.esa.int/documents/247904/1877131/DI-MPC-PB-0241-3-10_Sentinel-1IPFAuxiliaryProductSpecification.pdf/ae025687-c3e3-6ab0-de8d-d9cf58657431?t=1669115416469`

    Parameters
    ----------
    path_aux_cal: str 

    mode: str
        "IW" for example.

    pols : list 
        ["VV","VH"] for example;


    Returns
    ----------
    dict
        return a dict for the given (mode+pols).
        this dictionnary contains a dict with offboresight angle values and associated gains values
    """
    with open(path_aux_cal, 'rb') as file:
        xml = file.read()
        
    root_aux = objectify.fromstring(xml)
    dict_gains = {}

    for calibrationParams in root_aux.calibrationParamsList.getchildren():
        swath = calibrationParams.swath
        polarisation = calibrationParams.polarisation

        if (mode in swath.text and polarisation in pols):
            dict_temp = {}

            increment = calibrationParams.elevationAntennaPattern.elevationAngleIncrement
            valuesIQ = np.array([float(
                e) for e in calibrationParams.elevationAntennaPattern['values'].text.split(' ')])
            gain = np.sqrt(valuesIQ[::2]**2+valuesIQ[1::2]**2)

            count = gain.size
            ang = np.linspace(-((count - 1)/2) * increment,
                              ((count - 1)/2) * increment, count)

            dict_temp["offboresightAngle"] = ang
            dict_temp["gain"] = gain
            dict_gains[swath+"_"+polarisation] = dict_temp

    return dict_gains
