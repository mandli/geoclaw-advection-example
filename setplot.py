import os
import datetime

import numpy as np

import matplotlib
# # Markers and line widths
# matplotlib.rcParams['lines.linewidth'] = 2.0
# matplotlib.rcParams['lines.markersize'] = 6
# matplotlib.rcParams['lines.markersize'] = 8

# # Font Sizes
# matplotlib.rcParams['font.size'] = 16
# matplotlib.rcParams['axes.labelsize'] = 16
# matplotlib.rcParams['legend.fontsize'] = 12
# matplotlib.rcParams['xtick.labelsize'] = 16
# matplotlib.rcParams['ytick.labelsize'] = 16

# DPI of output images
article = False
if article:
    matplotlib.rcParams['savefig.dpi'] = 300
    matplotlib.rcParams['axes.titlesize'] = 'x-large'
    figsize_mult = 2
    add_colorbar = True
else:
    matplotlib.rcParams['savefig.dpi'] = 100
    figsize_mult = 1
    add_colorbar = True
import matplotlib.pyplot as plt

import clawpack.visclaw.colormaps as colormap
import clawpack.visclaw.gaugetools as gaugetools
import clawpack.clawutil.data as clawutil
import clawpack.amrclaw.data as amrclaw
import clawpack.geoclaw.data as geodata
import clawpack.geoclaw.surge.plot as surgeplot
import clawpack.geoclaw.topotools as topotools
import csv
from clawpack.geoclaw.util import fetch_noaa_tide_data

# surface_cmap = colormap.white_red
# surface_cmap = colormap.make_colormap({0:[0.0,1.0,1.0], 0.5:'#ffdddd', 1.0:"r"})
# surface_cmap = colormap.make_colormap({0:[0.0, 1.0, 1.0], 1.0:[1.0,0.0,0.0]})
# surface_cmap = colormap.make_colormap({0:'w', 1.0:'r'})

def setplot(plotdata=None):
    """"""

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    # clear any old figures,axes,items data
    plotdata.clearfigures()
    plotdata.format = 'ascii'

    # Load data from output
    clawdata = clawutil.ClawInputData(2)
    clawdata.read(os.path.join(plotdata.outdir, 'claw.data'))
    physics = geodata.GeoClawData()
    physics.read(os.path.join(plotdata.outdir, 'geoclaw.data'))
    surge_data = geodata.SurgeData()
    surge_data.read(os.path.join(plotdata.outdir, 'surge.data'))
    friction_data = geodata.FrictionData()
    friction_data.read(os.path.join(plotdata.outdir, 'friction.data'))

    # Load storm track
    track = surgeplot.track_data(os.path.join(plotdata.outdir, 'fort.track'))

    # Set afteraxes function
    def surge_afteraxes(cd):
        # 11:00 was 10:60 instead!
        surgeplot.surge_afteraxes(cd, track, plot_track=True, style="rX",
                                             kwargs={"markersize": 6})
        surgeplot.days_figure_title(cd, new_time=True)

    # Color limits
    surface_limits = [-1.0, 1.0]
    speed_limits = [0.0, 2.0]
    wind_limits = [0, 64]
    pressure_limits = [935, 1013]
    friction_bounds = [0.01, 0.04]

    def friction_after_axes(cd):
        plt.title(r"Manning's $n$ Coefficient")

    # ==========================================================================
    #   Plot specifications
    # ==========================================================================
    regions = {'Full Domain': {"xlimits": [clawdata.lower[0], clawdata.upper[0]],
                               "ylimits": [clawdata.lower[1], clawdata.upper[1]]}}

    # Surface Figure
    plotfigure = plotdata.new_plotfigure(name="Surface")
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Surface"
    plotaxes.xlimits = regions['Full Domain']['xlimits']
    plotaxes.ylimits = regions['Full Domain']['ylimits']
    plotaxes.scaled = True
    plotaxes.afteraxes = surge_afteraxes
    surgeplot.add_surface_elevation(plotaxes, bounds=surface_limits)
    surgeplot.add_land(plotaxes, bounds=[0.0, 20.0])
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0] * 10
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

    # Speed Figure
    plotfigure = plotdata.new_plotfigure(name="Currents")
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Currents"
    plotaxes.xlimits = regions['Full Domain']['xlimits']
    plotaxes.ylimits = regions['Full Domain']['ylimits']
    plotaxes.scaled = True
    plotaxes.afteraxes = surge_afteraxes
    surgeplot.add_speed(plotaxes, bounds=speed_limits)
    surgeplot.add_land(plotaxes, bounds=[0.0, 20.0])
    plotaxes.plotitem_dict['speed'].amr_patchedges_show = [0] * 10
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

    # Advective field
    plotfigure = plotdata.new_plotfigure(name="Advected Quantity")
    plotfigure.show = (clawdata.num_eqn > 3)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = r"$\psi$"
    plotaxes.xlimits = regions['Full Domain']['xlimits']
    plotaxes.ylimits = regions['Full Domain']['ylimits']
    plotaxes.scaled = True
    plotaxes.afteraxes = surge_afteraxes
    plotitem = plotaxes.new_plotitem(name='psi', plot_type='2d_pcolor')
    plotitem.plot_var = 4
    plotitem.pcolor_cmin = 0
    plotitem.pcolor_cmax = 1
    plotitem.pcolor_cmap = plt.get_cmap('cividis')
    plotitem.add_colorbar = True
    plotitem.colorbar_label = "Concentration"
    plotitem.amr_celledges_show = [0] * 10
    plotitem.amr_patchedges_show = [0] * 10
    surgeplot.add_land(plotaxes, bounds=[0.0, 20.0])
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

    # Friction field
    plotfigure = plotdata.new_plotfigure(name='Friction')
    plotfigure.show = friction_data.variable_friction and False
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = regions['Full Domain']['xlimits']
    plotaxes.ylimits = regions['Full Domain']['ylimits']
    plotaxes.title = "Manning's N Coefficient"
    plotaxes.afteraxes = friction_after_axes
    plotaxes.scaled = True
    surgeplot.add_friction(plotaxes, bounds=friction_bounds, shrink=0.9)
    plotaxes.plotitem_dict['friction'].amr_patchedges_show = [0] * 10
    plotaxes.plotitem_dict['friction'].colorbar_label = "$n$"

    # Pressure field
    plotfigure = plotdata.new_plotfigure(name='Pressure')
    plotfigure.show = surge_data.pressure_forcing and True
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = regions['Full Domain']['xlimits']
    plotaxes.ylimits = regions['Full Domain']['ylimits']
    plotaxes.title = "Pressure Field"
    plotaxes.afteraxes = surge_afteraxes
    plotaxes.scaled = True
    surgeplot.add_pressure(plotaxes, bounds=pressure_limits)
    surgeplot.add_land(plotaxes, bounds=[0.0, 20.0])
    plotaxes.plotitem_dict['pressure'].amr_patchedges_show = [0] * 10
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

    # Wind field
    plotfigure = plotdata.new_plotfigure(name='Wind Speed')
    plotfigure.show = surge_data.wind_forcing and True
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = regions['Full Domain']['xlimits']
    plotaxes.ylimits = regions['Full Domain']['ylimits']
    plotaxes.title = "Wind Field"
    plotaxes.afteraxes = surge_afteraxes
    plotaxes.scaled = True
    surgeplot.add_wind(plotaxes, bounds=wind_limits)
    surgeplot.add_land(plotaxes, bounds=[0.0, 20.0])
    plotaxes.plotitem_dict['wind'].amr_patchedges_show = [0] * 10
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    def get_actual_water_levels(station_id):
        # Fetch water levels and tide predictions for given station
        date_time, water_level, tide = fetch_noaa_tide_data(station_id,
                begin_date, end_date)

        # Calculate times relative to landfall
        days_rel_landfall = (date_time - landfall_time) / np.timedelta64(1,'s')  / 86400
        
        # Subtract tide predictions from measured water levels
        water_level -= tide
        
        # Find the mean values
        # Data imported every 6 minutes (i.e. 360 seconds)
        num_data_pts = (end_date - begin_date).total_seconds() / 360
        mean_value = np.sum(water_level) / num_data_pts
        # water_level -= mean_value
        print(f"{station_id}: {mean_value}")

        return days_rel_landfall, water_level

    def plot_observed(cd):
        station_id, station_name = stations[cd.gaugeno-1]
        days_rel_landfall, actual_level = get_actual_water_levels(station_id)

        ax = plt.gca()
        ax.plot(days_rel_landfall, actual_level, 'g')

    plotfigure = plotdata.new_plotfigure(name='Gauge Surfaces', figno=300,
                                         type='each_gauge')
    plotfigure.show = True
    plotfigure.clf_each_gauge = True

    stations = [ ('8518750', 'The Battery, NY'),
                 ('8516945', 'Kings Point, NY'),
                 ('8510560', 'Montauk, NY'),
                 ('8467150', 'Bridgeport, CT'),
                 ('8465705', 'New Haven, CT'),
                 ('8461490','New London, CT')]

    landfall_time = np.datetime64("2012-10-29T23:30")
    begin_date = datetime.datetime(2012, 10, 28, 23, 30)
    end_date = datetime.datetime(2012, 10, 30, 23, 30)
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.time_scale = 1 / (24 * 60**2)
    plotaxes.xlimits = [-1, 1]
    plotaxes.ylimits = [-1, 4]
    plotaxes.title = "Surface"
    plotaxes.ylabel = "Surface (m)"
    plotaxes.time_label = "Days relative to landfall"
    plotaxes.afteraxes = plot_observed
    plotaxes.grid = True

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    
    plotitem.plot_var = surgeplot.gauge_surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = surgeplot.gauge_dry_regions
    plotitem.kwargs = {"color":'lightcoral', "linewidth":5}

    #  Gauge Location Plot
    def gauge_location_afteraxes(cd):
        plt.subplots_adjust(left=0.12, bottom=0.06, right=0.97, top=0.97)
        surge_afteraxes(cd)
        gaugetools.plot_gauge_locations(cd.plotdata, gaugenos='all',
                                        format_string='ko', add_labels=True)

    plotfigure = plotdata.new_plotfigure(name="Gauge Locations")
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Gauge Locations'
    plotaxes.scaled = True
    plotaxes.xlimits = [-74.5, -71.6]
    plotaxes.ylimits = [40.5, 41.5]
    plotaxes.afteraxes = gauge_location_afteraxes
    surgeplot.add_surface_elevation(plotaxes, bounds=surface_limits)
    surgeplot.add_land(plotaxes, bounds=[0.0, 20.0])
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0] * 10
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

    # -----------------------------------------
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'   # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.latex = False                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # parallel plotting

    return plotdata