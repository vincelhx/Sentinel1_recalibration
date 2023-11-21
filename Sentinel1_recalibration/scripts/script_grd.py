import time

from Sentinel1_recalibration.utils import get_memory_usage, copy_tree, ignore_files
from Sentinel1_recalibration.S1_Recalibration import *

def recalibrate_grd(input_file, aux_version_config, config_file = None, resolution=None, save_netcdf=False, overwrite_netcdf=False, save_tiff=True, overwrite_tiff=False, fast = False):
    s1_recalibrer = S1_Recalibration(input_file, aux_version_config, config_file, resolution)
    
    if "sigma0_raw" not in INTEREST_VAR:
        logging.error(f"sigma0_raw must be in INTEREST_VAR={INTEREST_VAR}. Yet, we only use sigma0 to retrieve DN. Returning")
        return

    logging.info(f"save_netcdf is {save_netcdf} & save_tiff is {save_tiff}.")
    if save_netcdf:
        if resolution == None:
            s1_recalibrer.output_netcdf = os.path.join(
                s1_recalibrer.config["OUTPUTDIR"], "netcdf", aux_version_config, s1_recalibrer.SAFE, "output" +"fullres"+ ".nc")            
        else:
            s1_recalibrer.output_netcdf = os.path.join(
                s1_recalibrer.config["OUTPUTDIR"], "netcdf", aux_version_config, s1_recalibrer.SAFE, "output" +resolution+ ".nc")

        if fast:
            s1_recalibrer.output_netcdf = s1_recalibrer.output_netcdf.replace(".nc","_reduce.nc")
            
        if (overwrite_netcdf is False and os.path.exists(s1_recalibrer.output_netcdf)):
            logging.warn(
                f"save_netcdf -> netcdf file already exists and overwrite_netcdf is False.")
            save_netcdf = False

    if save_tiff:
        s1_recalibrer.outputdir_safe = os.path.join(
            s1_recalibrer.config["OUTPUTDIR"], "L1_recalibrated", aux_version_config ,s1_recalibrer.SAFE.replace(".SAFE", "_recal.SAFE"))
        
        
        s1_recalibrer.outputdir_measurement = os.path.join(
            s1_recalibrer.outputdir_safe, "measurement")

        measurement_vv = os.path.basename(
            glob.glob(os.path.join(s1_recalibrer.L1_path, "measurement/", "*-vv-*.tiff"))[0])
        measurement_vh = os.path.basename(
            glob.glob(os.path.join(s1_recalibrer.L1_path, "measurement/", "*-vh-*.tiff"))[0])

        s1_recalibrer.outputfile_vv = os.path.join(
            s1_recalibrer.outputdir_measurement, measurement_vv)
        s1_recalibrer.outputfile_vh = os.path.join(
            s1_recalibrer.outputdir_measurement, measurement_vh)

        if (overwrite_tiff is False and os.path.exists(s1_recalibrer.outputfile_vv) and os.path.exists(s1_recalibrer.outputfile_vh)):
            logging.warn(
                f"save_tiff -> vh & vv tiff files already exist and overwrite_tiff is False.")
            save_tiff = False

    if save_netcdf is False and save_tiff is False:
        logging.warn(
            "save_netcdf and save_tiff are now both False. Nothing to do. Returning.")
        return
    else:
        logging.info(f"save_netcdf is now {save_netcdf} & save_tiff is now {save_tiff}.")

    s1_recalibrer.create_interp_dict_geap()  # GRD only

    for pol in s1_recalibrer.polarizations:
        logging.info(f"starting pol={pol}")

        annot_slice = s1_recalibrer.get_annot_path(pol, 0)
        swath_bounds = s1_recalibrer.get_swath_bounds(annot_slice)

        offboresightAngle = s1_recalibrer.get_offboresight_angle(annot_slice)

        # swath_tab = np.zeros(dt.dataset.longitude.shape)#.astype(np.ubyte)
        swath_tab = np.full_like(
            s1_recalibrer.dataset.longitude, np.nan, dtype=float)
        for key in swath_bounds:
            swath_bound = swath_bounds[key]
            zloc = np.where((s1_recalibrer.dataset.line >= swath_bound['y_min']) & (s1_recalibrer.dataset.line <= swath_bound['y_max']) & (
                s1_recalibrer.dataset.sample >= swath_bound['x_min']) & (s1_recalibrer.dataset.sample <= swath_bound['x_max']))
            swath_tab[zloc] = int(swath_bound['swath'][-1])

        # Obtention de toutes les valeurs uniques dans swath_tab sans les nan et en int
        unique_slices = np.unique(swath_tab[~np.isnan(swath_tab)]).astype(int)

        #  Init Geap
        array_geap_old = np.full_like(offboresightAngle, np.nan, dtype=float)
        array_geap_new = np.full_like(offboresightAngle, np.nan, dtype=float)

        # Init Gproc
        array_gproc_old = np.full_like(offboresightAngle, np.nan, dtype=float)
        array_gproc_new = np.full_like(offboresightAngle, np.nan, dtype=float)

        for slice_value in unique_slices:
            
            logging.info(f"starting pol={pol}-slice={slice_value}")

            current_indices = np.where(swath_tab == slice_value)
            interp_geap_old = s1_recalibrer.dict_interp_geap_old.get(
                f'{s1_recalibrer.mode}{slice_value}_{pol}')
            array_geap_old[current_indices] = interp_geap_old(
                offboresightAngle.values[current_indices])

            interp_geap_new = s1_recalibrer.dict_interp_geap_new.get(
                f'{s1_recalibrer.mode}{slice_value}_{pol}')
            array_geap_new[current_indices] = interp_geap_new(
                offboresightAngle.values[current_indices])

            # SLC key_gproc = mode + "_" + product_type + "__1_" + mode + str(slice_value)
            key_gproc = s1_recalibrer.mode + "_" + \
                s1_recalibrer.product + "_1_" + \
                s1_recalibrer.mode + str(slice_value)
            infos_gproc_old = s1_recalibrer.dict_gproc_old[key_gproc]
            infos_gproc_old = {"HH": infos_gproc_old[0],
                               "HV": infos_gproc_old[1],
                               "VV": infos_gproc_old[2],
                               "VH": infos_gproc_old[3]}

            array_gproc_old[current_indices] = infos_gproc_old[pol]

            infos_gproc_new = s1_recalibrer.dict_gproc_new[key_gproc]
            infos_gproc_new = {"HH": infos_gproc_new[0],
                               "HV": infos_gproc_new[1],
                               "VV": infos_gproc_new[2],
                               "VH": infos_gproc_new[3]}

            array_gproc_new[current_indices] = infos_gproc_new[pol]

        Geap_old = xr.DataArray(
            array_geap_old, coords=offboresightAngle.coords)
        Geap_new = xr.DataArray(
            array_geap_new, coords=offboresightAngle.coords)

        Gproc_old = xr.DataArray(
            array_gproc_old, coords=offboresightAngle.coords)
        Gproc_new = xr.DataArray(
            array_gproc_new, coords=offboresightAngle.coords)

        s1_recalibrer.dataset["Geap_old_"+pol] = Geap_old
        s1_recalibrer.dataset["Gproc_old_"+pol] = Gproc_old
        s1_recalibrer.dataset["Geap_new_"+pol] = Geap_new
        s1_recalibrer.dataset["Gproc_new_"+pol] = Gproc_new

    for var in INTEREST_VAR:
        var_db = var+"_dB"
        s1_recalibrer.dataset[var_db] = 10*np.log10(s1_recalibrer.dataset[var])

        for pol in s1_recalibrer.polarizations:
            s1_recalibrer.dataset[var_db+"__"+pol+"__corrected"] = s1_recalibrer.dataset[var_db].sel(pol=pol) + 10*np.log10(s1_recalibrer.dataset["Geap_old_"+pol]) - 10*np.log10(s1_recalibrer.dataset["Geap_new_"+pol]) -\
                2*10*np.log10(s1_recalibrer.dataset["Gproc_old_"+pol]) + 2*10*np.log10(
                    s1_recalibrer.dataset["Gproc_new_"+pol])

            s1_recalibrer.dataset[var_db+"__"+pol +
                                  "__corrected"].expand_dims(pol=[pol])

        s1_recalibrer.dataset[var_db+"__corrected"] = xr.concat(
            [s1_recalibrer.dataset[var_db+"__"+'VV'+"__corrected"], s1_recalibrer.dataset[var_db+"__"+'VH'+"__corrected"]], dim='pol')

        s1_recalibrer.dataset[var+"__corrected"] = 10**(
            s1_recalibrer.dataset[var_db+"__corrected"]/10)

    if save_netcdf:
        s1_recalibrer.save_netcdf(fast)

    if save_tiff:

        if "sigma0_raw" in INTEREST_VAR:
            s1_recalibrer.s1dt._dataset["sigma0_raw"] = s1_recalibrer.dataset["sigma0_raw"+"__corrected"]

            dn_new = s1_recalibrer.s1dt.reverse_calibration_lut("sigma0_raw")
            s1_recalibrer.save_tiff(dn_new)
            copy_tree(s1_recalibrer.L1_path, s1_recalibrer.outputdir_safe, ignore=ignore_files)

        else:
            logging.error("sigma0_raw not in INTEREST_VAR. Returnin")
            return 
        
    s1_recalibrer.close_dss()


def processor_grd():
    import argparse
    import os
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description='Perform inversion from S1(L1-GRD) SAFE, L1-RCM, L1-RS2 ; using xsar/xsarsea tools')
    parser.add_argument('--input_file', help='input file path', required=True)
    
    parser.add_argument('--aux_version_config',
                        help='version name of the wanted new AUX conf files.',
                        required=False, default="v_IPF_36")
    
    parser.add_argument('--resolution',
                    help='resolution in km. only working if save_tiff is false',
                    required=False, default=None)
        
    parser.add_argument('--fast', action='store_true', default=False,
                        help='save in netcdf and only keep some vars',required=False)

    
    parser.add_argument('--verbose', action='store_true', default=False)

    parser.add_argument('--overwrite_netcdf', action='store_true', default=False,
                        help='overwrite existing .nc files [default is False]', required=False)

    parser.add_argument('--overwrite_tiff', action='store_true', default=False,
                        help='overwrite existing .tiff files [default is False]', required=False)

    parser.add_argument('--save_netcdf', action='store_true', default=False,
                        help='save as netcdf format [default is False]', required=False)

    parser.add_argument('--save_tiff', action='store_true', default=False,
                        help='save as tiff format [default is False]', required=False)
    
    parser.add_argument('--config_file', default=None
                        help='choose the path of the config file ; if None : default config', required=False)
    
    args = parser.parse_args()
    fmt = '%(asctime)s %(levelname)s %(filename)s(%(lineno)d) %(message)s'

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format=fmt,
                            datefmt='%d/%m/%Y %H:%M:%S', force=True)
    else:
        logging.basicConfig(level=logging.INFO, format=fmt,
                            datefmt='%d/%m/%Y %H:%M:%S', force=False)

    t0 = time.time()
    input_file = args.input_file.rstrip('/')
    
    logging.info('input file: %s', input_file)

        
    if args.resolution != None: 
        save_tiff = False
        logging.info('save_tiff set to False : we save_tiff at full resolution only')
        
    
    recalibrate_grd(input_file, args.aux_version_config, config_file = None, resolution = args.resolution, save_netcdf=args.save_netcdf,
                    overwrite_netcdf=args.overwrite_netcdf, save_tiff=args.save_tiff, overwrite_tiff=args.overwrite_tiff, fast=args.fast)

    logging.info('current memory usage: %s ', get_memory_usage(var='current'))
    logging.info('done in %1.3f min', (time.time() - t0) / 60.)


if __name__ == "__main__":
    processor_grd()
