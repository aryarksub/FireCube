import grpc
import os
import pyregence.fuel_wx_ign_pb2 as fuel_wx_ign_pb2
import pyregence.fuel_wx_ign_pb2_grpc as fuel_wx_ign_pb2_grpc
import requests
import tarfile

import util.general_util as gen_util
import util.processing_util as proc_util
  
def get_cloudfire_channel():
    if "CLOUDFIRE_SERVER" in os.environ:
        cloudfire_server = os.environ['CLOUDFIRE_SERVER']
    else:
        cloudfire_server ='worldgen.cloudfire.io'
    return cloudfire_server + ':50052'

def download_pyregence_data(
        out_dir, out_file, center, buffer=(60,60,60,60), wx_start_time=None, wx_num_hours=24, redo=False
):
    output_path = os.path.join(out_dir, out_file)
    if not redo and os.path.exists(output_path):
        return

    assert len(center) == 2 and len(buffer) == 4

    cloudfire_channel = get_cloudfire_channel()
    center_lat, center_lon = center
    west_buffer, east_buffer, south_buffer, north_buffer = [1000*buf for buf in buffer]

    with grpc.insecure_channel(cloudfire_channel) as channel:
        stub = fuel_wx_ign_pb2_grpc.FuelWxIgnStub(channel)
        response = stub.GetDomainData(fuel_wx_ign_pb2.Request( name = out_file ,
                                                               center_lat = center_lat,
                                                               center_lon = center_lon,
                                                               west_buffer = west_buffer,
                                                               east_buffer = east_buffer,
                                                               south_buffer = south_buffer,
                                                               north_buffer = north_buffer,
                                                               do_fuel = True,
                                                               fuel_source = 'landfire',
                                                               fuel_version = '2.4.0',
                                                               do_wx = True,
                                                               wx_type = 'historical',
                                                               wx_start_time = wx_start_time.strftime ("%Y-%m-%d %H:%M"),
                                                               wx_num_hours = wx_num_hours,
                                                               do_ignition = False,
                                                               point_ignition = True,
                                                               ignition_lat = -9999,
                                                               ignition_lon = -9999,
                                                               polygon_ignition = False,
                                                               active_fire_timestamp = None,
                                                               already_burned_timestamp = None,
                                                               ignition_radius = 300,
                                                               outdir = out_dir ) )

        # Stream the download to avoid loading the entire file into memory
        with requests.get(response.fileloc, stream=True) as r:
            r.raise_for_status()
            with open(output_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=1024*1024):
                    if chunk:  # filter out keep-alive chunks
                        f.write(chunk)

def extract_tif_from_pyr_tar(fid, tar_filepath):
    with tarfile.open(tar_filepath, 'r') as tar:
        for member in tar.getmembers():
            if member.name.endswith('tif'):
                var_name = os.path.splitext(os.path.basename(member.name))[0]
                output_path = gen_util.get_temp_data_video_filename(
                    fid, var_name, dir_type=gen_util.dir_data,
                    data_source=gen_util.subdir_pyr, var_type=gen_util.subdir_type_original
                )
                with tar.extractfile(member) as tif_file:
                    with open(output_path, 'wb') as out_file:
                        out_file.write(tif_file.read())

def driver_pyregence(fid, fire_center, fire_start, fire_hours, plot_types=[]):
    full_out_path = gen_util.get_pyr_tar_filename(fid)
    out_filename = os.path.basename(full_out_path)
    out_dir = os.path.dirname(full_out_path)
    
    download_pyregence_data(out_dir, out_filename, fire_center, (90,90,90,90), fire_start, fire_hours, redo=False)
    extract_tif_from_pyr_tar(fid, full_out_path)
    data_vars = gen_util.get_tif_vars_in_dir(
        os.path.join(gen_util.dir_temp, gen_util.dir_data, fid, gen_util.subdir_pyr, gen_util.subdir_type_original)
    )
    
    for var in data_vars:
        pyr_var_fnames = [
            gen_util.get_temp_data_video_filename(
                fid, var, dir_type=gen_util.dir_data,
                data_source=gen_util.subdir_pyr, var_type=vtype
            )
            for vtype in gen_util.var_types
        ]
        pyr_plot_fnames = [
            gen_util.get_temp_data_video_filename(
                fid, var, dir_type=gen_util.dir_videos,
                data_source=gen_util.subdir_pyr, var_type=vtype
            )
            for vtype in gen_util.var_types
        ]

        # Change CRS of original data to EPSG:5070
        proc_util.change_tif_crs(pyr_var_fnames[0], pyr_var_fnames[1], 'EPSG:5070')
        # Resample CRS-converted data to resolution defined in resample_tif (closest multiple of 30)
        proc_util.resample_tif(pyr_var_fnames[1], pyr_var_fnames[2], target_res=None)

        for plot_type in set(plot_types):
            if plot_type in gen_util.var_types:
                index = gen_util.var_types.index(plot_type)
                gen_util.create_animation_plot_from_tif(
                    in_tif=pyr_var_fnames[index],
                    out_file=pyr_plot_fnames[index],
                    start_time=fire_start,
                    mask=True,
                    ignore_small_neg=True
                )
