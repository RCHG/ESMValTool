"""

Diagnostics to estimate aerosols burden analysis.

Author: Ramiro Checa-Garcia (LSCE-IPSL)
        rcheca@lsce.ipsl.fr

Method:
    - It estimates global mean values and create time series
      with monthly and yearly time resolution.

Variables:
    - emidust, emisoa, emiss etc to estimate emission flux
    - drydust, drysoa, etc to estimate dry deposition flux

Outputs:
    - Single values or time series.

"""

import os
import numpy as np
import matplotlib
# use this everytime you import matplotlib
# modules; some machines dont have graphical interface (X)
matplotlib.use('Agg')  # noqa

import iris
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.preprocessor._area_pp import area_average

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
import logging
logger = logging.getLogger(os.path.basename(__file__))

## Specific modules needed by this diagnostic ----------------------------

import pandas as pd
from tabulate import tabulate
from pprint import pformat
#import xarray as xr
#daxr = xr.DataArray.from_iris(newcube)


def _get_my_files(cfg):
    """Put files in dicts of datasets and return them."""
    files_dict = {}
    for filename, attributes in cfg['input_data'].items():
        base_file = os.path.basename(filename)
        dataset = base_file.split('_')[1]
        files_dict[dataset] = {}
        files_dict[dataset]['file'] = filename
        if 'fx_files' in attributes:
            for fx_var in attributes['fx_files']:
                files_dict[dataset][fx_var] = attributes['fx_files'][fx_var]

    return files_dict

def _get_my_infos(cfg):
    """Put files in dicts of datasets and return them."""
    info_dict = {}
    for filename, attributes in cfg['input_data'].items():
        base_file = os.path.basename(filename)
        dataset = base_file.split('_')[1]
        info_dict[dataset] = {}
        info_dict[dataset]['units'] = attributes['units']
        info_dict[dataset]['short_name'] = attributes['short_name']
        if 'fx_files' in attributes:
            for fx_var in attributes['fx_files']:
                info_dict[dataset][fx_var] = attributes['fx_files'][fx_var]

    return info_dict


def main(cfg):
    """Compute the global specie emissions for each input dataset."""

    '''
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(input_data, 'standard_name', sort='dataset')
    logger.info( "Example of how to group and sort input data by standard_name:"
                 "\n%s", pformat(grouped_input_data))

    # Example of how to loop over variables/datasets in alphabetical order
    for standard_name in grouped_input_data:
        logger.info("Processing variable %s", standard_name)
        for attributes in grouped_input_data[standard_name]:
            logger.info("Processing dataset %s", attributes['dataset'])

            filename = attributes['filename']
            logger.info("Loading %s", filename)

            name = os.path.splitext(os.path.basename(filename))[0] + '_mean'
            logger.info("Name %s", name)
    '''

    global_emisions(cfg)

    return

def global_emisions(cfg):
    my_files_dict = _get_my_files(cfg)
    my_infos_dict = _get_my_infos(cfg)

    input_data = cfg['input_data'].values()

    print(my_infos_dict)
    grouped_input_data = group_metadata(
        input_data, 'standard_name', sort='dataset')
    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))

    # Example of how to loop over variables/datasets in alphabetical order
    for standard_name in grouped_input_data:
        logger.info("Processing variable %s", standard_name)
        for attributes in grouped_input_data[standard_name]:
            logger.info("Processing dataset %s", attributes['dataset'])
            filename = attributes['filename']
            logger.info("Loading %s", filename)
            fname = os.path.splitext(os.path.basename(filename))[0] + '_mean'
            fname = fname.replace(attributes['dataset'],'COMPARISON')
            logger.info("Name %s", fname)
            varname = attributes['short_name']
            logger.info("Units %s", attributes['units'])
    # Iterates through preprocessed model data to create multi-datasets tables

    gemissions_m = {}
    for key, value in my_files_dict.items():
        cube = iris.load_cube(value['file'])
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

        # Creates a new temporal cube with the global mean area values
        cube_area = iris.analysis.cartography.area_weights(cube)
        newcube = cube.collapsed(['longitude', 'latitude'], iris.analysis.SUM,
                                 weights=cube_area)

        # IMPORTANT --------------------------------------------------------- 
        # here we could provide a time series and seasons and sampled monthly

        try:
            tbounds = newcube.coord('time').bounds
            wbounds = [y[1]-y[0] for y in tbounds]
        except TypeError:
            newcube.coord('time').guess_bounds()
            tbounds = newcube.coord('time').bounds
            print(tbounds)
            wbounds = [y[1]-y[0] for y in tbounds]

        # we assume here time units of hours but variable in sec-1 ----------
        if my_infos_dict[key]['units']=='kg m-2 s-1':
            stime = np.array(wbounds)*24*60*60
            fdata = np.array([ a*b/1.e9  for a, b in zip(stime,newcube.data)])
            time = newcube.coord('time')
            atime = time.units.num2date(time.points)
            gemissions_m[key]=fdata
            gemissions_m['Date']=atime
        else:
            exit()

    globalemi_mon = pd.DataFrame(gemissions_m)
    globalemi_mon = globalemi_mon.set_index('Date')
    globalemi_yrs = globalemi_mon.resample("Y").sum()

    if cfg.get('table'):
        path_table_mon = os.path.join(
                    cfg['table_dir'],
                    fname + '_mon.' +  cfg['table_output_file_type'],
                    )
        path_table_yrs = os.path.join(
                    cfg['table_dir'],
                    fname + '_yrs.' +  cfg['table_output_file_type'],
                    )
        path_table_des = os.path.join(
                    cfg['table_dir'],
                    fname + '_des.' +  cfg['table_output_file_type'],
                    )

        if not os.path.exists(cfg['table_dir']):
            os.makedirs(cfg['table_dir'])

        if cfg['output_file_type']!='md':
            month_tb = tabulate(globalemi_mon, headers='keys', tablefmt='psql')
            years_tb = tabulate(globalemi_yrs, headers='keys', tablefmt='psql')
            month_de = tabulate(globalemi_mon.describe(), headers='keys', tablefmt='psql')
            years_de = tabulate(globalemi_yrs.describe(), headers='keys', tablefmt='psql')

        if cfg['output_file_type']!='tex':
            month_tb = tabulate(globalemi_mon, headers='keys', tablefmt='latex')
            years_tb = tabulate(globalemi_yrs, headers='keys', tablefmt='latex')
            month_de = tabulate(globalemi_mon.describe(), headers='keys', tablefmt='latex')
            years_de = tabulate(globalemi_yrs.describe(), headers='keys', tablefmt='latex')

        print(path_table_mon, path_table_yrs, path_table_des)
        if cfg['table']['monthly']==True:
            with open(path_table_mon, 'w') as tablef:
                tablef.write(month_tb)
        if cfg['table']['yearly']==True:
            with open(path_table_yrs, 'w') as tablef:
                tablef.write(years_tb)
        if cfg['table']['summary']==True:
            with open(path_table_des, 'w') as tablef:
                tablef.write(month_de)
                tablef.write(years_de)

    if cfg.get('plot'):
        path_plot_mon = os.path.join(
                    cfg['plot_dir'],
                    fname + '_mon.' +  cfg['output_file_type'],
                    )
        path_plot_yrs = os.path.join(
                    cfg['plot_dir'],
                    fname + '_yrs.' +  cfg['output_file_type'],
                    )
        if cfg['plot']['monthly']==True:
            globalemi_mon.plot(figsize=(12,4), legend=True)
            plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.90), facecolor=None, edgecolor=None, frameon=False)
            plt.title('Comparison of global monthly '+varname)
            plt.ylabel(varname + ' [Tg month-1]')
            plt.subplots_adjust(right=0.7)
            plt.savefig(path_plot_mon)
        if cfg['plot']['yearly']==True:
            globalemi_yrs.plot(figsize=(12,4), legend=True)
            plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.90), facecolor=None, edgecolor=None, frameon=False)
            plt.title('Comparison of global yearly '+varname)
            plt.ylabel(varname +' [Tg yr-1]')
            plt.subplots_adjust(right=0.7)
            plt.savefig(path_plot_yrs)

        #plt.plot(newcube.coord('time'), fdata)
        #plt.show()
        # Dado un intervalo de tiempo largo es deseable que saque:
        # median and mean sobre todos los anos
        # median and mean sobre todos los meses
        # extremos -> dar year y mes
        # pdfs
        # time series
        # mean seasonal values
        # anomalies
        #
        # here we have to ascertain the number of seconds per month and then use this
        # to create a yearly value, or a monthly value.
    return


def plot_time_series(cfg):
    """
    Example of personal diagnostic function.

    Arguments:
        run - dictionary of data files

    Returns:
        string; makes some time-series plots

    """
    # local path for e.g. plots: user input
    root_dir = '/group_workspaces/jasmin2/cmip6_prep/'  # edit as per need
    out_path = 'esmvaltool_users/valeriu/'   # edit as per need
    local_path = os.path.join(root_dir, out_path)

    # get the files (simple case, one-variable only)
    my_files_dict = _get_my_files(cfg)

    # iterate through preprocessed model data
    for key, value in my_files_dict.items():
        cube = iris.load_cube(value['file'])
        area_avg_cube = area_average(cube, 'latitude', 'longitude')
        plt.plot(area_avg_cube.data[:, 0], label=key)
        plt.xlabel('Time (months)')
        plt.ylabel(cube.standard_name)
        plt.title('Time series at ground level')
        plt.tight_layout()
        plt.grid()
        plt.legend()
        png_name = 'Time_series_' + key + '.png'
        plt.savefig(os.path.join(local_path, png_name))
        plt.close()

    return 'I made some plots!'


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)

