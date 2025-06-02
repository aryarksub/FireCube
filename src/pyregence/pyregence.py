import grpc
import os
import pyregence.fuel_wx_ign_pb2 as fuel_wx_ign_pb2
import pyregence.fuel_wx_ign_pb2_grpc as fuel_wx_ign_pb2_grpc
import requests

import util.general_util as gen_util
  
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

def driver_pyregence(fid, fire_center, fire_start, fire_hours):
    full_out_path = gen_util.get_pyr_tar_filename(fid)
    out_filename = os.path.basename(full_out_path)
    out_dir = os.path.dirname(full_out_path)
    
    download_pyregence_data(out_dir, out_filename, fire_center, (90,90,90,90), fire_start, fire_hours, redo=True)
