from lxml import etree, objectify
import numpy as np
import xarray as xr

def get_bounds(annot_file_path):
    """
    register the swath merging information
    """
    from shapely.geometry import Polygon

    tree = etree.parse(annot_file_path)
    swath_bounds_dict = dict()

    mode = tree.xpath('//adsHeader/swath')[0].text
    if mode == 'EW':
        beam_list = ['EW1', 'EW2', 'EW3', 'EW4', 'EW5']
    elif mode == 'IW':
        beam_list = ['IW1', 'IW2', 'IW3']
    else:
        raise NotImplementedError()

    for beam in beam_list:
        pattern = f'//swathMergeList/swathMerge/swathBoundsList/swathBounds[../../swath/text()="{beam}"]'
        for sb_number, swath_bounds in enumerate(tree.xpath(pattern)):
            x_min = int(swath_bounds.find('firstRangeSample').text)
            x_max = int(swath_bounds.find('lastRangeSample').text)
            y_min = int(swath_bounds.find('firstAzimuthLine').text)
            y_max = int(swath_bounds.find('lastAzimuthLine').text)

            # Define the corners
            corners = [(x_min, y_min), (x_min, y_max),
                       (x_max, y_max), (x_max, y_min)]
            # Create the polygon
            poly = Polygon(corners)

            # Store in data structure
            swath_bounds_dict[(beam, sb_number)] = {
                'swath': beam,
                'x_min': x_min,
                'x_max': x_max,
                'y_min': y_min,
                'y_max': y_max,
                'geometry': poly}

    return swath_bounds_dict


# ** FCTS to gets OffborsightAngle **/

def get_dict_roll_slc(annot):
    """
    find rollAngles(azimuthTime)

    Parameters
    ----------
    annot: str 
        path of the annotation file

    Returns
    ----------
    dict
        return a dict key,value = azimuthTime, rollAngles
    """
    # Parse the XML file
    with open(annot, 'rb') as file:
        xml = file.read()

    root = objectify.fromstring(xml)
    dict_roll = {}
    for antennaPattern in root.antennaPattern.antennaPatternList.getchildren():
        dict_roll[antennaPattern.azimuthTime.text] = float(antennaPattern.roll)
    return dict_roll


def getRollAngle(dates, interp):
    """
    find roll angles associated with dates for a given function 'interp'

    Parameters
    ----------
    dates: np.ndarray
        1d datetime64[ns]

    interp: scipy.interpolate._interpolate.interp1d
        relation between time and roll angles obtained from get_dict_roll_slc values

    Returns
    ----------
    np.ndarray
        roll angle values
    """
    dates_sec = dates.astype('int64') / 1e9
    return interp(dates_sec)


def getOffboresightAngle(annot, aztime, elevation):
    """
    get offboresight angle values
    1/ get roll values from annot file,
    2/ convert the dates values in seconds to create a function interp(time[s]) = roll_angle[s]
    3/ use this function with the aztime dates from the L1. with getRollAngle
    4/ return offboresight_angle = elevation - roll_angles

    Parameters
    ----------
    annot: str 
        path of the annotation file
    aztim: np.ndarray 
        azimuth time values, 1d datetime64[ns]
    elevation: np.ndarray 
        elevation angle values, 1d float64

    Returnsf
    ----------
    np.ndarray
        offboresight angle values, 1d float64
    """
    from scipy.interpolate import interp1d
    from datetime import datetime
    import time
    dict_roll = get_dict_roll_slc(annot)

    timestamps = list(dict_roll.keys())

    def convert_to_seconds(t): return time.mktime(datetime.fromisoformat(
        t).timetuple()) + datetime.fromisoformat(t).microsecond / 1e6
    seconds = np.array(list(map(convert_to_seconds, timestamps)))

    # relation between time and roll angles
    interp = interp1d(seconds, np.array(list(dict_roll.values())),
                      kind='linear', fill_value='extrapolate')
    # flatten the array and convert it to numpy array
    aztime_np = aztime.values.ravel()
    # apply getRollAngle to the flattened array
    roll_angles = getRollAngle(aztime_np, interp)
    # reshape the result and convert it back to xr.DataArray
    roll_angles = xr.DataArray(roll_angles.reshape(
        aztime.shape), coords=aztime.coords)
    offboresight_angle = elevation - roll_angles

    return offboresight_angle
