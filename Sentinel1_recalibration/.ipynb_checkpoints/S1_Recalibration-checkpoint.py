import os
import logging
import glob
import numpy as np
import xsar
import xarray as xr
import yaml

from Sentinel1_recalibration.aux_paths import get_aux_paths
from Sentinel1_recalibration.aux_pp1_tools import get_all_gproc
from Sentinel1_recalibration.aux_cal_tools import get_geap_gains
from Sentinel1_recalibration.annot_tools import get_bounds, getOffboresightAngle

def load_config():
    with open("./../config.yaml", 'r') as stream:
        config = yaml.safe_load(stream)
    return config

config = load_config()
AUX_CAL = config["AUX_CAL"]
AUX_PP1 = config["AUX_PP1"]
PATH_DATA = config["PATH_DATA"]
INTEREST_VAR = ["sigma0_raw"]
OUTPUTDIR = config["OUTPUTDIR"]


class S1_Recalibration:
    def __init__(self, L1_path, resolution=None):

        self.config = config

        self.L1_path = L1_path
        self.SAFE = os.path.basename(L1_path)

        self.s1dt = xsar.Sentinel1Dataset(L1_path, resolution=resolution)

        self.product_type = self.s1dt.sar_meta.manifest_attrs['product_type']
        self.product = self.s1dt.sar_meta.product

        self.mode = self.s1dt.sar_meta.manifest_attrs['swath_type']
        self.polarizations = self.s1dt.sar_meta.manifest_attrs["polarizations"].tolist(
        )

        self.resolution = resolution
        self.dataset = self.s1dt.dataset

        self.get_aux_paths()
        self.get_new_AUX_CAL_path()
        self.get_new_AUX_PP1_path()

        self.dict_geap_old = self.get_geap_dict(self.PATH_AUX_CAL_OLD)
        self.dict_geap_new = self.get_geap_dict(self.PATH_AUX_CAL_NEW)

        self.dict_gproc_old = self.get_gproc_dict(self.PATH_AUX_PP1_OLD)
        self.dict_gproc_new = self.get_gproc_dict(self.PATH_AUX_PP1_NEW)

    def get_aux_paths(self):
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
        aux_paths = get_aux_paths(self.L1_path)
        self.PATH_AUX_CAL_OLD = PATH_DATA+"/AUX_CAL/" + \
            aux_paths['AUX_CAL'] + "/data/" + \
            aux_paths['AUX_CAL'][0:3].lower()+"-aux-cal.xml"
        self.PATH_AUX_PP1_OLD = PATH_DATA+"/AUX_PP1/" + \
            aux_paths['AUX_PP1'] + "/data/" + \
            aux_paths['AUX_PP1'][0:3].lower()+"-aux-pp1.xml"
        # PATH_AUX_INS_OLD = path_original_data+"AUX_INS/" + aux_paths['AUX_INS'] +"/data/" + aux_paths['AUX_INS'][0:3].lower()+"-aux-ins.xml"

    def get_annot_path(self, pol: str, slice_number: int):
        if self.product_type == "GRD":
            return glob.glob(
                self.L1_path+"/annotation/*"+self.mode.lower()+"*"+pol.lower()+"*")[0]
        else:
            return glob.glob(self.L1_path+"/annotation/*"+self.mode.lower()+str(slice_number).lower()+"*"+pol.lower()+"*")[0]

    def get_new_AUX_CAL_path(self):
        self.PATH_AUX_CAL_NEW = PATH_DATA+"/AUX_CAL/"+self.config["AUX_CAL"]+"/data/" + \
            self.SAFE[0:3].lower()+"-aux-cal.xml"

    def get_new_AUX_PP1_path(self):
        self.PATH_AUX_PP1_NEW = PATH_DATA+"/AUX_PP1/"+self.config["AUX_PP1"]+"/data/" + \
            self.SAFE[0:3].lower()+"-aux-pp1.xml"

    def get_geap_dict(self, path_aux_cal):
        return get_geap_gains(path_aux_cal, self.mode, self.polarizations)

    def get_gproc_dict(self, path_aux_pp1):
        return get_all_gproc(path_aux_pp1)

    def get_swath_bounds(self, annot_file_path):
        return get_bounds(annot_file_path)

    def get_offboresight_angle(self, annot_file_path):
        return getOffboresightAngle(annot_file_path, self.dataset["azimuth_time"], self.dataset["elevation"])

    def create_interp_dict_geap(self):
        from scipy.interpolate import interp1d

        ret = []
        for _dict in [self.dict_geap_old, self.dict_geap_new]:
            # New dictionary to store the results
            curr_dict_interp_geap = {}
            # Loop through each key in the initial dictionary
            for key in _dict:
                # Get the information for the current key
                infos_geap = _dict[key]
                # Create the linear interpolation
                interp_geap = interp1d(
                    infos_geap["offboresightAngle"], infos_geap['gain'], kind='linear')
                # Add the interpolation to the new dictionary with the same key
                curr_dict_interp_geap[key] = interp_geap
            ret.append(curr_dict_interp_geap)

        # Assign the interpolated dictionaries to their respective attributes
        self.dict_interp_geap_old, self.dict_interp_geap_new = ret

    def save_netcdf(self, overwrite):
        if (overwrite or not os.path.exists(self.output_netcdf)):

            for var in INTEREST_VAR:
                var_db = var+"_dB"
                dataset = self.dataset.copy()
                dataset = dataset.drop(var_db+"__"+'VV'+"__corrected")
                dataset = dataset.drop(var_db+"__"+'VH'+"__corrected")

            dataset = dataset.drop_vars(["spatial_ref", "digital_number"])

            for var in ['footprint', 'multidataset', 'rawDataStartTime', 'specialHandlingRequired']:
                if var in dataset.attrs:
                    dataset.attrs[var] = str(dataset.attrs[var])
                if "approx_transform" in dataset.attrs:
                    del dataset.attrs["approx_transform"]

            os.makedirs(os.path.dirname(self.output_netcdf), exist_ok=True)
            dataset.to_netcdf(self.output_netcdf, mode="w")
            dataset.close()

        else:
            logging.info(
                f"save_netcdf -> File {self.output_netcdf} already exists and overwrite is False.")

    def save_tiff(self, dn_new, overwrite):

        data_array_VV = dn_new.sel(pol="VV").rio.set_spatial_dims(
            y_dim="line", x_dim="sample")
        data_array_VH = dn_new.sel(pol="VH").rio.set_spatial_dims(
            y_dim="line", x_dim="sample")

        # For VV
        data_array_VV.rio.to_raster(self.outputfile_vv)
        # For VH
        data_array_VH.rio.to_raster(self.outputfile_vh)

    def close_dss(self):
        self.dataset.close()
        self.s1dt.dataset.close()
        logging.info(
            f"dataset & s1dt.dataset closed")
