# Mitchell Chandler
# mitchell_chandler@tws.org
# Last updated: 16-March-2026
##################################################

#import Pkg

using NetCDF
using CairoMakie
using ColorSchemes
using Statistics
using HypothesisTests
using DataFrames
using GLM 

cd("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/_Data/scripts/results")

#Each grid cell is 10-km x 10-km = 100-km2

## READ IN DATA ##
# GAP
gap = ncread("ogdp50_gap.nc","grid_code");
gx = ncread("ogdp50_gap.nc","x");
gy = ncread("ogdp50_gap.nc","y");

# non-GAP
nongap = ncread("ogdp50_nongap.nc","grid_code");
nx = ncread("ogdp50_nongap.nc","x");
ny = ncread("ogdp50_nongap.nc","y");

# Plot all points
f = Figure();
ax = Axis(f[1, 1],aspect = DataAspect())
sc = scatter!(gx,gy,color=gap,colorrange=(0,1),colormap=:YlOrBr,markersize=3)
sc = scatter!(nx,ny,color=nongap,colorrange=(0,1),colormap=:YlOrBr,markersize=3)
Colorbar(f[1,2],sc,ticks=([0:0.1:1;]),height=300,width=20,tellheight=false)
hidedecorations!(ax)
display(f)

## COMPARE GAP V. NON-GAP LANDS ##
total_land = [gap; nongap];

# Find areas 
total_land_area = length(total_land)*100 #km2
gap_land_area = length(gap)*100; #km2

pct_gap_land = length(gap)/length(total_land)*100

# Mean and median
mm = DataFrame(
    Class=["GAP","non-GAP"],
    Mean=[mean(gap), mean(nongap)],
    Median=[median(gap), median(nongap)]
)

# Kolmogorov-Smirnov Test
ApproximateTwoSampleKSTest(gap,nongap)


## AREA OF HIGH OGDP ##
#terciles
low_ogdp_thresh = quantile(total_land,1/3);
high_ogdp_thresh = quantile(total_land,2/3)

high_ogdp_area = length(total_land[total_land.>=high_ogdp_thresh])*100; #km2

HO_ng_area = length(nongap[nongap.>=high_ogdp_thresh])*100; #km2
HO_gap_area = length(gap[gap.>=high_ogdp_thresh])*100; #km2

pct_HO_gap_land = HO_gap_area/high_ogdp_area*100
pct_HO_ng_land = HO_ng_area/high_ogdp_area*100


## PERCENT OF PROTECTED LANDS AT EACH PROBABILITY THRESHOLD ##
prob_threshold = [0.95:-0.05:0;];
all_pct = NaN*prob_threshold;
for i in eachindex(prob_threshold)
    all_pct[i] = length(findall(gap.>=prob_threshold[i]))./length(findall(total_land.>=prob_threshold[i]))*100;
end

## FIGURE ##
function plot_linear(x,y)
    #plot data
    lines!(x,y,color=:black,linewidth=4)

    #linear fit
    df = DataFrame(; x,y);
    model = lm(@formula(y ~ x),df);
    slope_p = coeftable(model).cols[4][2];

    #plot trend if significant
    if slope_p<0.05
        lines!(x,predict(model),color=:red,linewidth=2,linestyle=(:dash,:dense))
    end

    return coeftable(model)
end

f = Figure();

# PDF
ax = Axis(f[1,1],xlabel="OGDP",ylabel="Probability Density",
    xgridvisible=false,ygridvisible=false,
    xminorticks=IntervalsBetween(2),xminorticksvisible=true,xminorgridvisible=false)
#plot
density!(ax,gap,label="Protected",color=(Makie.to_colormap(:BrBG_3)[3],0.7),strokecolor=Makie.to_colormap(:BrBG_3)[3],strokewidth=3)
density!(ax,nongap,label="Unprotected",color=(Makie.to_colormap(:BrBG_3)[1],0.7),strokecolor=Makie.to_colormap(:BrBG_3)[1],strokewidth=3)
#axis:
xlims!(0,1)
ylims!(0,4)
ax.xticks=0.1:0.2:1.1
#legend:
axislegend(framevisible=false,position=:rt,backgroundcolor=:white,alpha=0,orientation=:vertical,patchsize=(30,15),rowgap=5)
#text:
text!(0,1,text="(a)",align=(:left,:top),offset=(2,-4),space=:relative)

# Trend
ax = Axis(f[1,2],
    xlabel="OGDP threshold",ylabel="Percent protected",
    xgridvisible=false,ygridvisible=false,
    xminorticks=IntervalsBetween(2),xminorticksvisible=true,xminorgridvisible=false,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true,yminorgridvisible=false)
#plot
vlines!(low_ogdp_thresh,color=:grey30,linestyle=(:dot,:dense),linewidth=2)
vlines!(high_ogdp_thresh,color=:grey30,linestyle=(:dot,:dense),linewidth=2)
plot_linear(prob_threshold,all_pct)
#axis:
xlims!(0,1)
ylims!(10,50)
ax.xticks=0.1:0.2:1.1
ax.yticks=0:10:50
ax.ytickformat=y -> string.(Int64.(y)).*"%"
#text:
text!(0,1,text="(b)",align=(:left,:top),offset=(2,-4),space=:relative)

display(f)
#save("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/Figures/ogdp_availability-2026.02.12.pdf",f)

