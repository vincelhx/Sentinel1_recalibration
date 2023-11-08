from lxml import objectify
import os


def get_aux_paths(L1_path):
    """
    Find the paths of the auxiliary files associated with the L1 file. 

    https://github.com/umr-lops/xsar/issues/174 could make it easier 

    Parameters
    ----------
    L1_path: str 

    Returns
    ----------
    dict
        return a dict of auxiliary file paths 

    """

    aux_paths = {}

    path_manifest = L1_path + "/manifest.safe"

    # Parse the XML file
    with open(path_manifest, 'rb') as file:
        xml = file.read()
    root_manifest = objectify.fromstring(xml)

    # Define the namespaces used in the XML
    namespaces = {
        'safe': 'http://www.esa.int/safe/sentinel-1.0',
    }

    # Find all resource elements
    resource_elements = root_manifest.findall(
        './/safe:resource', namespaces=namespaces)

    # Get name and role attributes of each resource and select the AUX ones
    for resource in resource_elements:
        role = resource.get('role')
        if (role in ["AUX_PP1", "AUX_CAL", "AUX_INS"]):
            name = resource.get('name')
            aux_paths[role] = os.path.basename(name)
    return aux_paths


""" TODO WHEN https://github.com/umr-lops/xsar/issues/174 solved
def get_auxPaths(L1_path):
    "return aux_cal_path and aux_pp1_path used" 
    #aux_cal_path = os.path.basename(meta.manifest_attrs['aux_cal'])
    
    aux_paths = {
    "AUX_PP1" = os.path.basename(meta.manifest_attrs['aux_cal'])
    "AUX_CAL" = os.path.basename(meta.manifest_attrs['aux_pp1'])
    "AUX_INS" = os.path.basename(meta.manifest_attrs['aux_ins'])
    }
"""
