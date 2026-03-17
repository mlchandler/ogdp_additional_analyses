# Mitchell Chandler
# mitchell_chandler@tws.org
# Last updated: 13-February-2026
############################

using CairoMakie
using CSV
using DataFrames
using Statistics

cd("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/_Data/scripts/eia")

## Read in data ##
world_oil = CSV.read("eia_world_crude_oil_lease_condensate_prod_2026.02.13.csv",DataFrame,transpose=true); #Mb/d = thousand barrels per day
world_gas = CSV.read("eia_world_dry_natural_gas_prod_2026.02.13.csv",DataFrame,transpose=true); #bcf = billion cubic feet


#Convert oil from thousand barrels per day -> barrels
world_oil[:,:bbl] = world_oil[:,:World]*365.25*1E3

#Convert gas from billion cubic feet -> trillion cubic feet
world_gas[:,:tcf] = world_gas[:,:World]*1E9/1E12;


#Plot
f = Figure(fontsize=16);
#oil
ax1 = Axis(f[1,1],ylabel="Billion barrels",
    ygridvisible=false,xgridvisible=false,
    title="(a) Global oil production",titlealign=:right,titlefont=:regular)
#data - barrels
lines!(world_oil[:,:Year],world_oil[:,:bbl]/1E9,linewidth=5,color=:black)
#axis:
xlims!(1970,2025)
ax1.xticks=1970:10:2025;
ylims!(15,35)

#gas
ax2 = Axis(f[2,1],ylabel="Trillion cubic feet",
    ygridvisible=false,xgridvisible=false,
    title="(b) Global gas production",titlealign=:right,titlefont=:regular)
#data - trillion cubic feet
lines!(world_gas[:,:Year],world_gas[:,:tcf],linewidth=5,color=:black)
#axis:
xlims!(1970,2025)
ax2.xticks=1970:10:2025;
ylims!(40,160)
ax2.yticks=40:40:160;

ax2.yticklabelspace = maximum(tight_yticklabel_spacing!,[ax1,ax2]);
ax1.yticklabelspace = maximum(tight_yticklabel_spacing!,[ax1,ax2]);

display(f)
#save("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/Figures/global_og_production-2026.02.13.pdf",f)


#Calculations
#2024 production
oil2024 = world_oil[world_oil[:,:Year].==2024,:bbl][1];
oil2024/1e9
gas2024 = world_gas[world_gas[:,:Year].==2024,:tcf][1]

#percent increase from 1980 to 2024
oil1980 = world_oil[world_oil[:,:Year].==1980,:bbl][1];
gas1980 = world_gas[world_gas[:,:Year].==1980,:tcf][1];

oil_pct_increase = (oil2024-oil1980)/oil1980*100
gas_pct_increase = (gas2024-gas1980)/gas1980*100

