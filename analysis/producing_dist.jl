# Mitchell Chandler
# mitchell_chandler@tws.org
# Last updated: 13-April-2026
##################################################

#import Pkg

using CairoMakie
using Statistics
using NetCDF

cd("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/_Data/scripts/results")

## READ IN DATA ##  
#OGDP at well locations
wells_ogdp = ncread("../50pct_witheld/ogdp50_composite_prediction_points.nc","grid_code")

#OGDP at background locations
bground_ogdp = ncread("../50pct_witheld/ogdp50_composite_background_points.nc","grid_code");

ogdp = [wells_ogdp; bground_ogdp];

## PERCENT OF PRODUCING CELLS CONSIDERED HIGH OGDP ##
#Terciles
low_ogdp_thresh = quantile(ogdp,1/3);
high_ogdp_thresh = quantile(ogdp,2/3);

length(findall(wells_ogdp.>=high_ogdp_thresh)) / length(wells_ogdp) * 100

## PLOT PDF ##
f = Figure();
ax = Axis(f[1,1],xlabel="OGDP",ylabel="Probability Density",
    xgridvisible=false,ygridvisible=false,
    xminorticks=IntervalsBetween(2),xminorticksvisible=true,xminorgridvisible=false)
#plot
density!(ax,wells_ogdp,color=(:black,0.3),strokecolor=:black,strokewidth=1)
vlines!(high_ogdp_thresh,color=:grey30,linestyle=(:dot,:dense),linewidth=2)
#axis:
xlims!(0,1)
ylims!(0,2.5)
ax.xticks=0.1:0.2:1.1
display(f)
