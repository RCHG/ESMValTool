"""
Autoassess Consservation - energy budget.

Module with routines to calculate energy conservation on various sub-models.
Presently, it only includes the atmospheric energy budget
It will include other sub-models in the future.
"""

import matplotlib.pyplot as plt
import numpy as np
import cf_units
import iris
import iris.analysis.calculus as icalc

from esmvaltool.diag_scripts.autoassess.loaddata import load_run_ss
from esmvaltool.preprocessor._area_pp import area_average_general as a_avg


def atmos_energy_budget(run):
    """
    Compute budget.

    Function to check atmospheric energy budgets and conservation.
    It uses section 30 diagnostics.
    It uses Mean diagnostics (fluxes) and
    initial and final fields for the meaning period

    Program originally written in IDL for New Dynamics  by  R A Stratton.
    10/01/2015 porting of program to IRIS by Jose Rodriguez
    """
    metrics = dict()

    # constants required for calculations
    # use AuxCoords in order to maintain units
    l_c = iris.coords.AuxCoord(
        2.501e6, long_name='lantent_heat_of_condensation', units='J kg-1')
    l_f = iris.coords.AuxCoord(
        0.334e6, long_name='lantent_heat_of_fusion', units='J kg-1')
    secs_per_season = iris.coords.AuxCoord(
        90. * 86400.0, long_name='seconds_per_season', units='s')

    # missing-data indicator:
    # TODO use NaN for missing data
    mdi = -10000.0

    # Make sure we pick up last instantaneous field
    endyear = run['to_annual']
    endyear = endyear.replace(year=endyear.year + 1)

    try:
        # Switch from stash to standard_name
        # column integral cvT per unit area
        # cvt_f = load_run_ss(run, 'instantaneous',
        #                     'm01s30i420', to_dt=endyear)
        # column integral gr per unit area
        # gr_f  = load_run_ss(run, 'instantaneous',
        #                     'm01s30i421', to_dt=endyear)
        # TOTAL KE PER UA WITH W  RHO GRID
        # ke_f  = load_run_ss(run, 'instantaneous',
        #                     'm01s30i402', to_dt=endyear)
        # energy correction
        # en_cor = load_run_ss(run, 'seasonal', 'm01s30i419')

        # Seasonal means of energy fluxes:

        # m01s01i207: incoming SW rad flux (TOA)
        # CMOR: rsdt, exact long name
        swin = load_run_ss(run, 'seasonal', 'toa_incoming_shortwave_flux')

        # m01s01i208: outgoing SW rad flux (TOA)
        # CMOR: rsut, exact name
        swout = load_run_ss(run, 'seasonal', 'toa_outgoing_shortwave_flux')

        # m01s01i201: net down surface SW flux
        # CMOR: rsds, surface_downwelling_shortwave_flux_in_air
        # orig diag name: surface_net_downward_shortwave_flux
        sw_flux = load_run_ss(run, 'seasonal',
                              'surface_downwelling_shortwave_flux_in_air')
        # m01s02i201: net down surface LW flux
        # CMOR: rlds, surface_downwelling_longwave_flux_in_air
        # orig diag name: surface_net_downward_longwave_flux
        lw_flux = load_run_ss(run, 'seasonal',
                              'surface_downwelling_longwave_flux_in_air')
        # m01s03i332: TOA outgoing LW rad
        # CMOR: rlut, exact name
        olr = load_run_ss(run, 'seasonal', 'toa_outgoing_longwave_flux')

        # m01s03i217: surface heat flux"
        # CMOR: hfss, exact name
        sh_flux = load_run_ss(run, 'seasonal',
                              'surface_upward_sensible_heat_flux')

        # m01s05i215: total snowfall rate
        # CMOR: prsn, exact name
        snow = load_run_ss(run, 'seasonal', 'snowfall_flux')

        # m01s05i216: total precipitation rate
        # CMOR: pr, exact name
        precip = load_run_ss(run, 'seasonal', 'precipitation_flux')

        # VPREDOI TODO
        # temporarily replacing missing variables with pr
        # cvt, gr, ke, en_cor missing correct var defs
        cvt_f = load_run_ss(run, 'monthly', 'precipitation_flux', lbproc=192)
        gr_f = cvt_f  # missing
        ke_f = cvt_f  # missing
        en_cor = cvt_f  # missing
        cvt_f, gr_f, ke_f = remove_forecast_period([cvt_f, gr_f, ke_f])

        # Set appropriate units for fields loaded above
        cvt_f.units = cf_units.Unit('J m-2')
        gr_f.units = cf_units.Unit('J m-2')
        ke_f.units = cf_units.Unit('J m-2')

        # Remove forecast periods
        swin, swout, sw_flux, lw_flux, olr, sh_flux, snow, precip, en_cor = \
            remove_forecast_period([swin, swout,
                                    sw_flux, lw_flux, olr, sh_flux,
                                    snow, precip, en_cor])

        # calculate global budgets

        # calculate global averages and budgets:
        # instantaneous fields
        cvtg = a_avg(cvt_f, weighted=True)
        grg = a_avg(gr_f, weighted=True)
        keg = a_avg(ke_f, weighted=True)
        en_tot = keg + cvtg + grg

        # Rate of change of instantaneous fields
        ch_en = icalc.cube_delta(en_tot, 'time') / secs_per_season
        ch_ke = icalc.cube_delta(keg, 'time') / secs_per_season
        ch_cvt = icalc.cube_delta(cvtg, 'time') / secs_per_season
        ch_gr = icalc.cube_delta(grg, 'time') / secs_per_season

        # Energy fluxes
        swing = a_avg(swin, weighted=True)
        swoutg = a_avg(swout, weighted=True)
        swg = a_avg(sw_flux, weighted=True)
        lwg = a_avg(lw_flux, weighted=True)
        olrg = a_avg(olr, weighted=True)
        shg = a_avg(sh_flux, weighted=True)
        snowg = a_avg(snow, weighted=True)
        precipg = a_avg(precip, weighted=True)
        en_corg = a_avg(en_cor, weighted=True)

        # energy flux into atmosphere = radTOA - SH +
        # Lc * precip + Lf * snowfall
        toa = swing - swoutg - olrg
        diab_heat = toa - swg - lwg + shg + precipg * l_c + snowg * l_f

        # Remove time bounds of cube with time
        # average data to allow arithmetics
        # with cube containing instantaneous time points
        diab_heat.coord('time').bounds = None

        # VPREDOI TODO use this instead when units match
        # this is due to using pr instead of actual diag variable
        # err_en = ch_en - diab_heat
        err_en = ch_en

        ecorrection = np.mean(en_corg.data)
        eerror = np.mean(err_en.data)

        # assign metric:
        metrics['atmospheric global energy error'] = np.abs(eerror)

        # Produce extra plot:
        expid = run['runid']
        fig = plt.figure(figsize=(8.27, 11.69))

        x_err = np.arange(err_en.data.size)
        stitl1_temp = '{0} mean energy-conservation error: {1:7.4f} W/m2'
        stitl1 = stitl1_temp.format(expid, eerror)
        stitl2 = 'Energy correction: {0:7.4f} W/m2'.format(
            ecorrection)

        plt.subplot(2, 1, 1)
        titl1 = 'Change in total energy over 3 months'
        plt.plot(
            x_err, ch_en.data, linewidth=2,
            color='black', label='E(end)-E(start)')
        # plt.plot(x_err, diab_heat.data, linewidth=2, color='red',
        #          label='E added to atmos from fluxes')
        plt.plot(
            x_err, ch_ke.data, linewidth=2,
            color='blue', label='Change in KE')
        plt.plot(
            x_err, ch_cvt.data, linewidth=2,
            color='green', label='Change in cvT')
        plt.plot(
            x_err, ch_gr.data, linewidth=2,
            color='purple', label='Change in gr')
        plt.xlim([0, err_en.data.size])
        plt.ylim([np.min(ch_en.data) * 1.5, np.max(ch_en.data) * 1.5])
        plt.title(titl1)
        plt.xlabel('No. of seasons from djf ' +
                   str(run['from_annual'].year + 1))
        plt.ylabel('W/m2')
        plt.axhline(0.0, linestyle=':', color='black')
        plt.legend(loc='lower center', fontsize='small', frameon=0)

        plt.subplot(2, 1, 2)
        titl1 = 'Error in energy conservation over 3 month periods ' + expid
        plt.plot(x_err, err_en.data, linewidth=2, color='black')
        plt.xlim([0, err_en.data.size])
        plt.ylim([np.min(err_en.data) * 1.5, np.max(err_en.data) * 1.5])
        plt.title(titl1)
        plt.xlabel('No. of seasons from djf ' +
                   str(run['from_annual'].year + 1))
        plt.ylabel('W/m2')
        plt.axhline(0.0, linestyle=':', color='black')

        plt.suptitle(stitl1 + '\n' + stitl2, fontsize='large')
        plt.savefig(expid + '_atmospheric_energy_budget.png')

    except iris.exceptions.ConstraintMismatchError as msg:
        print(msg)
        metrics['atmospheric global energy error'] = mdi

    return metrics


def _remove_forecast_period(cube):
    try:
        cube.remove_coord('forecast_period')
    except iris.exceptions.CoordinateNotFoundError as exc:
        print(exc)
    return cube


def remove_forecast_period(cubes):
    """Remove forecast period."""
    return [_remove_forecast_period(cube) for cube in cubes]