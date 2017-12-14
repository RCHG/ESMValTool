diag_script_info = True

;; Used by: perfmetrics_main.ncl

;; Plot type ("cycle" (time), "zonal" (plev, lat), "latlon" (lat, lon))
diag_script_info@plot_type = "cycle"

;; Time average ("opt" argument of time_operations.ncl)
diag_script_info@time_avg = "monthlyclim"

;; Selected region ("Global", "Tropics", "NH extratropics", "SH extratropics")
diag_script_info@region = "Global"

;; Draw plots
diag_script_info@draw_plots = True

;; Draw difference plots
diag_script_info@plot_diff = True

;; Plot standard deviation ("all", "none", "ref_model" or given model name)
diag_script_info@plot_stddev = "ref_model"

;; Plot legend in a separate file (for cycle)
diag_script_info@legend_outside = True

;; Plot style (for cycle)
diag_script_info@styleset = "CMIP5"

;; Calculate t-test in difference plot (for zonal and lat-lon)
diag_script_info@t_test = False

;; Confidence level for the t-test (for zonal and lat-lon)
diag_script_info@conf_level = 0.95