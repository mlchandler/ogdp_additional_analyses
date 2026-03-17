# Mitchell Chandler
# mitchell_chandler@tws.org
# Last updated: 06-January-2026
##################################################

#import Pkg

using NetCDF
using CairoMakie
using ColorSchemes

cd("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/_Data/scripts/50pct_witheld")

##################################################
## READ IN DATA ##

# Bias (Copeland minus OGDP) at all points
bias = ncread("diff_points50.nc","grid_code");
biasx = ncread("diff_points50.nc","x");
biasy = ncread("diff_points50.nc","y");
#plot all bias points
f = Figure();
ax = Axis(f[1, 1],aspect = DataAspect())
sc = scatter!(biasx,biasy,color=bias,colorrange=(-1,1),colormap=Reverse(:RdYlBu),markersize=3)
Colorbar(f[1,2],sc,ticks=([-1:0.1:1;]),height=300,width=20,tellheight=false)
hidedecorations!(ax)
display(f)

# Bias (Copeland minus OGDP) at well points
wellbias = ncread("diff_prediction_points50.nc","grid_code");

#Copeland probabilities
copeland_prob = ncread("overlapping_copeland_points50.nc","grid_code");
copelandx = ncread("overlapping_copeland_points50.nc","x");
copelandy = ncread("overlapping_copeland_points50.nc","y");
copeland_well_prob = ncread("overlapping_copeland_prediction_points50.nc","grid_code"); 
#plot all Copeland points
f = Figure();
ax = Axis(f[1, 1],aspect = DataAspect())
sc = scatter!(copelandx,copelandy,color=copeland_prob,colorrange=(0,1),colormap=:YlOrBr,markersize=3)
Colorbar(f[1,2],sc,ticks=([-1:0.1:1;]),height=300,width=20,tellheight=false)
hidedecorations!(ax)
display(f)

#OGDP probabilities
ogdp_prob = ncread("overlapping_ogdp_points50.nc","grid_code");
ogdp_well_prob = ncread("overlapping_ogdp_prediction_points50.nc","grid_code"); 

##################################################
## PLOT ##
f = Figure();

# Histogram of all points (coloured bins)
ax = Axis(f[1, 1],ylabel="Probability density",
    xgridvisible=false,ygridvisible=false,
    xminorticks=IntervalsBetween(5),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true)
hist!(ax,bias,bins=(-1:0.1:1),normalization=:pdf,color=Makie.to_colormap(reverse(resample(ColorSchemes.RdYlBu,20))),strokecolor=:grey10,strokewidth=0.1)
vlines!(0,color=:grey60,linestyle=:solid,linewidth=1)
#text:
text!(0,1,text="(b)",align=(:left,:top),offset=(3,-1),space=:relative)
#axis:
xlims!(-1,1)
ylims!(0,2)
ax.yticks=0:1:2;

# Histogram of well points (coloured bins)
ax = Axis(f[2, 1],ylabel="Probability density",xlabel="Bias",
    xgridvisible=false,ygridvisible=false,
    xticksmirrored=true,
    xminorticks=IntervalsBetween(5),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true)
hist!(ax,wellbias,bins=(-1:0.1:1),normalization=:pdf,color=Makie.to_colormap(reverse(resample(ColorSchemes.RdYlBu,20))),strokecolor=:grey20,strokewidth=0.1)
vlines!(0,color=:grey60,linestyle=:solid,linewidth=1)
#text:
text!(0,1,text="(c)",align=(:left,:top),offset=(3,-1),space=:relative)
#axis:
xlims!(-1,1)
ylims!(0,2)
ax.yticks=0:1:2;

# PDF of all points
ax = Axis(f[1, 2],
    xgridvisible=false,ygridvisible=false,
    xminorticks=IntervalsBetween(5),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true)
#OGDP:
density!(ax,ogdp_prob,label="OGDP",color=(Makie.to_colormap(:RdYlBu_3)[3],0.7),strokecolor=Makie.to_colormap(:RdYlBu_3)[3],strokewidth=2)
#Copeland:
density!(ax,copeland_prob,label="Copeland",color=(Makie.to_colormap(:RdYlBu_3)[1],0.7),strokecolor=Makie.to_colormap(:RdYlBu_3)[1],strokewidth=2)
#text:
text!(0,1,text="(d)",align=(:left,:top),offset=(3,-1),space=:relative)
#axis:
xlims!(0,1)
ylims!(0,6)
ax.yticks=0:1:6;
#legend:
axislegend(framevisible=false,position=:rt,backgroundcolor=:white,alpha=0,patchsize=(15,15),orientation=:vertical)

# PDF of well points
ax = Axis(f[2, 2],xlabel="Probability",
    xgridvisible=false,ygridvisible=false,
    xticksmirrored=true,
    xminorticks=IntervalsBetween(5),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true)
#OGDP:
density!(ax,ogdp_well_prob,label="OGDP",color=(Makie.to_colormap(:RdYlBu_3)[3],0.7),strokecolor=Makie.to_colormap(:RdYlBu_3)[3],strokewidth=2)
#Copeland:
density!(ax,copeland_well_prob,label="Copeland",color=(Makie.to_colormap(:RdYlBu_3)[1],0.7),strokecolor=Makie.to_colormap(:RdYlBu_3)[1],strokewidth=2)
#text:
text!(0,1,text="(e)",align=(:left,:top),offset=(3,-1),space=:relative)
#axis:
xlims!(0,1)
ylims!(0,4)
ax.yticks=0:1:4;
#legend:
axislegend(framevisible=false,position=:rt,backgroundcolor=:white,alpha=0,patchsize=(15,15),orientation=:vertical) #labelsize=10

display(f)
#save("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/Figures/subimages/copeland_comparison-julia-2026.01.06.svg",f) #note this svg needs to be optimised to be opened in ArcGIS Pro e.g. using SVGOMG

##################################################
# Find mode of a probability density function
function pdf_mode(data)
    k = Makie.KernelDensity.kde(data);
    k_mode = k.x[argmax(k.density)];
    return k_mode
end

pdf_mode(bias)
pdf_mode(ogdp_well_prob)
pdf_mode(copeland_well_prob)

