# Mitchell Chandler
# mitchell_chandler@tws.org
# Last updated: 23-January-2026
##################################################

#import Pkg

using CairoMakie
using CSV
using DataFrames
using Statistics
using GLM 

cd("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/_Data/scripts/results")

prob_threshold = [0.95:-0.05:0;];

##################################################
## READ IN DATA ##  
df = CSV.read("ogdp_slope_etopo_gap.csv",DataFrame);
ogdp = df[:,:grid_code];
slope = df[:,:degrees];
elevation = df[:,:metres];
protected = df[:,:protected];

##################################################
## COMPUTATIONS FOR ALL GRID CELLS GREATER THAN OR EQUAL TO A GIVEN OGDP THRESHOLD ##  
av_slope = NaN*prob_threshold;
av_elevation = NaN*prob_threshold;
frac_protected_o = NaN*prob_threshold;
for i in eachindex(prob_threshold)
    idx = findall(ogdp.>=prob_threshold[i]);
    av_slope[i] = mean(slope[idx]);
    av_elevation[i] = mean(elevation[idx]);
    frac_protected_o[i] = length(findall(protected[idx].==1))/length(idx);    
end

##################################################
## OGDP THRESHOLD AGAINST AVERAGE SLOPE ##
#Linear fit
model = lm(@formula(av_slope ~ prob_threshold),DataFrame(; prob_threshold,av_slope))

#Plot
f = Figure();
ax = Axis(f[1,1],
    xlabel="Probability threshold",ylabel="Average slope [degrees]",title="All",
    xminorticks=IntervalsBetween(2),xminorticksvisible=true,xminorgridvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true,yminorgridvisible=true)
lines!(prob_threshold,av_slope,color=:black,linewidth=4)
lines!(prob_threshold,predict(model),color=:red,linewidth=2,linestyle =:dash)
#axis:
xlims!(0,1)
ylims!(0,1)
ax.xticks=0:0.2:1
ax.yticks=0:0.2:1
display(f)

##################################################
## COMPUTATIONS FOR GRID CELLS WITHIN A SLOPE BIN ##  
slope_edges = collect(0:0.25:ceil(maximum(slope)));
slope_centres = (slope_edges[1:end-1]+slope_edges[2:end])/2;
av_ogdp = NaN*slope_centres;
frac_protected_s = NaN*slope_centres;
for i in eachindex(slope_centres)
    idx = findall(slope_edges[i] .<= slope .< slope_edges[i+1]);
    av_ogdp[i] = mean(ogdp[idx]);
    frac_protected_s[i] = length(findall(protected[idx].==1))/length(idx);
end

##################################################
## BINNED SLOPE AGAINST AVERAGE OGDP ##
# Remove NaNs
slope_centres2 = slope_centres[.!isnan.(av_ogdp)];
av_ogdp2 = av_ogdp[.!isnan.(av_ogdp)];

# Linear fit
model = lm(@formula(y ~ x),DataFrame(x=slope_centres2,y=av_ogdp2))

# Plot
f = Figure();
ax = Axis(f[1,1],
    xlabel="Slope [degrees]",ylabel="Average OGDP",title="All",
    xgridvisible=false,ygridvisible=false)
barplot!(slope_centres2,av_ogdp2,gap=0,color=:grey20,strokecolor=:gray60,strokewidth=2)
lines!(slope_edges,Float64.(predict(model,DataFrame(x=slope_edges))),color=:red,linewidth=3,linestyle =:dash)
#axis:
xlims!(minimum(slope_edges),maximum(slope_edges))
ylims!(0,0.6)
ax.xticks=0:1:8
ax.yticks=0:0.1:1
display(f)

##################################################
## BINNED SLOPE AGAINST FRACTION PROTECTED ##
# Remove NaNs
slope_centres2 = slope_centres[.!isnan.(frac_protected_s)];
frac_protected2 = frac_protected_s[.!isnan.(frac_protected_s)];

# Linear fit
model = lm(@formula(y ~ x),DataFrame(x=slope_centres2,y=frac_protected2))

# Plot
f = Figure();
ax = Axis(f[1,1],
    xlabel="Slope [degrees]",ylabel="Fraction protected",title="All",
    xgridvisible=false,ygridvisible=false)
barplot!(slope_centres2,frac_protected2,gap=0,color=:grey20,strokecolor=:gray60,strokewidth=2)
lines!(slope_edges,Float64.(predict(model,DataFrame(x=slope_edges))),color=:red,linewidth=3,linestyle =:dash)
#axis:
xlims!(minimum(slope_edges),maximum(slope_edges))
ylims!(0,1)
ax.xticks=0:1:8
ax.yticks=0:0.1:1
display(f)
