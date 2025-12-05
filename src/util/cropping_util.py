import rasterio
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import imageio_ffmpeg as imgffm
import matplotlib as mpl
from matplotlib.patches import Rectangle

ffmpeg_path = imgffm.get_ffmpeg_exe()
mpl.rcParams['animation.ffmpeg_path'] = ffmpeg_path

plt.switch_backend('agg')

from processing_util import pad_bounds_to_resolution_multiple
import general_util as gen_util

def determine_band_range(fperim_tif):
    """
    Determine the first and last time bands containing fire perimeters. Since there are 1-day buffers
    when computing bounding boxes and bands are hourly time-steps, the first and last 24 bands should
    not be counted towards the first/last time bands.

    Args:
        fperim_tif (str): Path to the fire perimeter TIF file.

    Returns:
        tuple: A tuple containing the first and last time bands with fire perimeters.
    """
    # Minimum time band = first band where there is a fire perimeter
    # Maximum time band = last band where there is a fire perimeter
    # Note: there are 1-day buffers when computing bounding boxes (min band > 24, max band <= src.count-24)
    with rasterio.open(fperim_tif) as src:
        min_band = src.count
        max_band = 0
        for band in range(25, src.count - 24 + 1):
            data = src.read(band)
            if np.any(data == 1):
                min_band = min(min_band, band)
                max_band = max(max_band, band)

    return min_band, max_band

def get_bounding_box_per_timestep(fperim_tif, coarsest_res=9000):
    """
    Compute spatial bounding boxes for each timestep with fire perimeters. The bounding boxes
    are padded so that the height and width are both multiples of the given coarsest resolution.
    The bounding boxes are stored as tuples of the form (minx, miny, maxx, maxy). 

    Args:
        fperim_tif (str): Path to the fire perimeter TIF file.
        coarsest_res (int, optional): The coarsest resolution (in meters) to which bounding boxes 
        should be padded. Defaults to 9000.

    Returns:
        dict: Dictionary mapping each time band to its corresponding bounding box tuple (minx, miny, maxx, maxy).
    """
    bboxes = dict()
    # Determine first/last time bands with fire perimeters
    min_band, max_band = determine_band_range(fperim_tif)
    
    with rasterio.open(fperim_tif) as src:
        transform = src.transform

        for band in range(min_band, max_band + 1):
            data = src.read(band)
            
            # For each band, determine which pixels represent burned area
            burned_area_mask = (data == 1)

            rows, cols = np.where(burned_area_mask)

            min_row, max_row = rows.min(), rows.max()
            min_col, max_col = cols.min(), cols.max()

            # convert pixel indices to bounding box in spatial coords
            # UL (upper-left) has max y since spatial grid is north-up/south-down
            minx, maxy = rasterio.transform.xy(transform, min_row, min_col, offset='ul')
            maxx, miny = rasterio.transform.xy(transform, max_row, max_col, offset='lr')

            # Pad bounding box to be multiples of coarsest resolution
            final_x1, final_y1, final_x2, final_y2 = pad_bounds_to_resolution_multiple(
                (minx, miny, maxx, maxy), coarsest_res=coarsest_res
            )

            bboxes[band] = (final_x1, final_y1, final_x2, final_y2)
    
    return bboxes

def restrict_tifs_to_bbox(in_tifs, out_tifs, bounds_by_band):
    """
    Restrict each of the input TIFs to a bounding box for each time band given in the bounds_by_band
    dictionary. The data inside the bounding box will be kept while all data outside the bounding
    box is replaced by np.nan (NaN) values. Single-band input TIFs are written to output as 
    multi-band TIFs with data taken from the single band but restricted to bounding boxes as
    appropriate. Each input TIF has a corresponding output TIF where the data will be saved. Note
    that due to the functionality for single-band input TIFs expanding to multi-band TIFs, the total 
    size of the output TIFs will likely be much larger than the total size of the input TIFs.

    Args:
        in_tifs (list): List of paths to the input TIF files. Each input TIF can be single-band or multi-band.
        out_tifs (list): List of paths to the output TIF files. Each output TIF corresponds to an input TIF.
        bounds_by_band (dict): Dictionary mapping each time band to its corresponding bounding box tuple 
         (minx, miny, maxx, maxy).
    """
    assert len(in_tifs) == len(out_tifs) and len(bounds_by_band) > 0
    for in_tif, out_tif in zip(in_tifs, out_tifs):
        with rasterio.open(in_tif) as src:
            # Copy input profile and update count/compression for output TIF
            profile = src.profile.copy()
            height, width = src.height, src.width
            profile.update(
                count=len(bounds_by_band),
                dtype='float32',
                nodata=np.nan,
                compress='lzw'
            )

            output_bands = []
            for band in bounds_by_band:
                # x/y min/max spatial bounds
                xmin, ymin, xmax, ymax = bounds_by_band[band]
                # Band index needs to be chosen appropriately for multi-band vs single-band TIFs
                band_index = band if src.count > 1 else 1
                data = src.read(band_index).astype('float32')

                # Convert x/y spatial bounds to pixel indices
                # Upper-left (min row/col) has max y since spatial grid is north-up/south-down
                row_min, col_min = src.index(xmin, ymax)
                # Lower-right (max row/col) has min y since spatial grid is north-up/south-down
                row_max, col_max = src.index(xmax, ymin)

                # Clip to raster extent to avoid out-of-bounds
                row_min = max(0, min(row_min, height-1))
                row_max = max(0, min(row_max, height-1))
                col_min = max(0, min(col_min, width-1))
                col_max = max(0, min(col_max, width-1))

                # Create mask for bounding box and set data outside bounding box to NaN
                mask = np.zeros_like(data, dtype=bool)
                mask[row_min:row_max+1, col_min:col_max+1] = True
                masked_data = data.copy()
                masked_data[~mask] = np.nan
                output_bands.append(masked_data)
    
        # Write output TIF with all bands
        with rasterio.open(out_tif, 'w', **profile) as dst:
            for i, band_data in enumerate(output_bands, start=1):
                dst.write(band_data, i)

def driver(fids=[], plot=False, redo=False): 
    """
    Driver function for TIF cropping and optional batch visualization.

    Args:
        fids (list, optional): List of fire IDs to process. If an empty list is given, then process
         all stored fires. Defaults to [].
        plot (bool, optional): True if animation plots should be created; False otherwise. Defaults to False.
        redo (bool, optional): True to reprocess all fires even if output exists; False to skip existing fires. 
         Defaults to False.
    """
    # Get all fids in the output/cubes and run cropping on each with possibility for animation
    fids = fids if len(fids) > 0 else gen_util.get_all_output_fids()

    for fid_index, fid in enumerate(fids, start=1):
        fperim_tif = os.path.join(gen_util.dir_output, gen_util.dir_cubes, fid, gen_util.subdir_firespread, 'fperim.tif')

        # Get bounding boxes per time band based on fire perimeter
        bboxes = get_bounding_box_per_timestep(fperim_tif)
        band_range = determine_band_range(fperim_tif)
        shifted_bboxes = {i : bboxes[i + band_range[0] - 1] for i in range(1, len(bboxes) + 1)}

        for batch in gen_util.data_batches:
            in_tifs = []
            out_tifs = []
            batch_exists = False

            if batch != gen_util.subdir_vis:
                # Original cube output directory
                batch_dir = os.path.join(
                    gen_util.dir_output, gen_util.dir_cubes, fid, batch
                )
                # Cropped cube output directory (create if it doesn't exist)
                out_dir = os.path.join(
                    gen_util.dir_output, gen_util.dir_crop_cubes, fid, batch
                )
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                else:
                    batch_exists = True
                # Visualization output directory (create if it doesn't exist)
                vis_dir = os.path.join(
                    gen_util.dir_output, gen_util.dir_crop_cubes, fid, gen_util.subdir_vis
                )
                if not os.path.exists(vis_dir):
                    os.makedirs(vis_dir)

                animation_fname = os.path.join(vis_dir, f'{batch}.mp4')
                
                for var_file in os.listdir(batch_dir):
                    if var_file.endswith('.tif'):
                        in_tifs.append(os.path.join(batch_dir, var_file))
                        out_tifs.append(os.path.join(out_dir, f'cropped_{var_file}'))

                # Perform cropping if redo is True or batch doesn't already exist
                if redo or not batch_exists:
                    restrict_tifs_to_bbox(
                        in_tifs=in_tifs,
                        out_tifs=out_tifs,
                        bounds_by_band=bboxes
                    )
                
                # Create animation if requested
                if plot:
                    animate_multiple_tifs_with_bboxes(
                        tif_paths=out_tifs,
                        bbox_dict=shifted_bboxes,
                        out_file=animation_fname
                    )
        
        if fid_index % 10 == 0:
            print(f'Processed cropping for {fid_index} fires')

def bbox_to_pixel_coords(ds, bbox):
    """
    Convert a spatial bounding box to pixel coordinates for a raster dataset.

    Args:
        ds (rasterio.io.DatasetReader): An open raster dataset.
        bbox (tuple): A bounding box in the form (minx, miny, maxx, maxy).

    Returns:
        tuple: Pixel coordinates in the form (x0, y0, width, height).
    """
    minx, miny, maxx, maxy = bbox

    # Convert world to pixel coordinates using raster transform
    row_min, col_min = ds.index(minx, maxy)  # upper-left
    row_max, col_max = ds.index(maxx, miny)  # lower-right

    x0 = col_min
    y0 = row_min
    width = col_max - col_min
    height = row_max - row_min

    return x0, y0, width, height

def animate_multiple_tifs_with_bboxes(
    tif_paths,
    bbox_dict,
    out_file,
    num_cols=4,
    figsize_scale=4,
    cmap="viridis",
    fps=6
):
    """
    Create an animation showing multiple TIF rasters with per-band bounding boxes. Each input TIF is displayed 
    in its own subplot. The animation advances through time by iterating over bands shared across all rasters. 
    For every band, the corresponding bounding box is drawn as an overlaid rectangle.

    Args:
        tif_paths (list): List of paths to the input TIF files to be animated.
        bbox_dict (dict): Dictionary mapping band indices to bounding boxes.
        out_file (str): Path to the output animation file.
        num_cols (int, optional): Number of subplot columns. Defaults to 4.
        figsize_scale (int, optional): Multiply figure width/height by this factor. Defaults to 4.
        cmap (str, optional): Colormap. Defaults to "viridis".
        fps (int, optional): Frames per second for the animation. Defaults to 6.
    """
    # Load datasets and determine global time steps
    datasets = [rasterio.open(p) for p in tif_paths]
    band_list = list(bbox_dict.keys())

    num_frames = len(band_list)

    # Compute subplot grid
    n = len(tif_paths)
    num_rows = int(np.ceil(n / num_cols))

    fig, axes = plt.subplots(
        num_rows, num_cols,
        figsize=(figsize_scale * num_cols, figsize_scale * num_rows)
    )
    axes = np.array(axes).reshape(-1)

    # Turn off unused axes
    for j in range(n, len(axes)):
        axes[j].axis("off")

    # Initial setup per subplot
    im_list = []
    rect_list = []

    for ax, ds, tif_path in zip(
        axes[:n], datasets, tif_paths
    ):
        first_band = band_list[0]
        arr0 = ds.read(first_band)

        # Compute bounding box rectangle for first band
        x0, y0, w, h = bbox_to_pixel_coords(ds, bbox_dict[first_band])

        # Show image
        im = ax.imshow(arr0, cmap=cmap)
        im_list.append(im)

        # Add rectangle for bounding box
        rect = Rectangle(
            (x0, y0),
            w, h,
            fill=False,
            edgecolor="red",
            linewidth=2
        )
        ax.add_patch(rect)
        rect_list.append(rect)

        # Title
        title = os.path.basename(tif_path).replace(".tif", "")
        ax.set_title(f"{title} - {first_band}")

    # Update function for animation
    def update(frame_idx):
        for plot_idx in range(n):
            ds = datasets[plot_idx]

            # If fewer bands than max_t, freeze last frame
            if frame_idx < len(band_list):
                band = band_list[frame_idx]
            else:
                band = band_list[-1]

            arr = ds.read(band)
            im_list[plot_idx].set_data(arr)
            im_list[plot_idx].set_clim(np.nanmin(arr), np.nanmax(arr))

            # Update bounding box
            x0, y0, w, h = bbox_to_pixel_coords(ds, bbox_dict[band])
            rect_list[plot_idx].set_xy((x0, y0))
            rect_list[plot_idx].set_width(w)
            rect_list[plot_idx].set_height(h)

            # Update title
            title = os.path.basename(tif_paths[plot_idx]).replace(".tif", "")
            axes[plot_idx].set_title(f"{title} - {band_list[frame_idx]}")

        return im_list + rect_list

    # Construct animation
    anim = FuncAnimation(
        fig,
        update,
        frames=num_frames,
        interval=1000 // fps,
        blit=False
    )

    # Save animation
    writer = FFMpegWriter(fps=fps)
    anim.save(out_file, writer=writer)
    plt.close(fig)


if __name__=='__main__':
    fid = 'TX3000510039120150811'

    driver([], plot=False, redo=False)