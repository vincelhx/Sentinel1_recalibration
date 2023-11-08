from lxml import objectify
import numpy as np


def get_all_gproc(path_aux_pp1):
    """
    Find gains Gproc associated with mode product and slice number from AUX_PP1.

    DOC : `https://sentinel.esa.int/documents/247904/1877131/DI-MPC-PB-0241-3-10_Sentinel-1IPFAuxiliaryProductSpecification.pdf/ae025687-c3e3-6ab0-de8d-d9cf58657431?t=1669115416469`

    Parameters
    ----------
    path_aux_pp1: str 

    Returns
    ----------
    dict
        return a dict of 4 linear gain values for each (mode+product_type+slice_number)
        in this order : Gproc_HH, Gproc_HV, Gproc_VV, Gproc_VH
    """

    # Parse the XML file
    with open(path_aux_pp1, 'rb') as file:
        xml = file.read()
    root_pp1 = objectify.fromstring(xml)
    dict_gains = {}
    for product in root_pp1.productList.getchildren():
        for swathParams in product.slcProcParams.swathParamsList.getchildren():
            gains = [float(g) for g in swathParams.gain.text.split(' ')]
            dict_gains[product.productId+"_"+swathParams.swath] = gains
    return dict_gains
