# Mitchell Chandler
# mitchell_chandler@tws.org
# Last updated: 03-February-2026
##################################################

#import Pkg

using CSV
using DataFrames

cd("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/_Data/scripts/results")

##################################################
## READ IN DATA ##  
df = CSV.read("bivariate_points.csv",DataFrame);
biv = df[:,:grid_code];
state = Vector(df[:,:STATE_ABBR]);

#for a given bivariate category, compute the percentage of that category in a given state 
#note that bivariate category is in the form [ogdp tercile][wcv tercile] e.g. 13 is low ogdp & high wcv
function pct_cat_state(biv_cat,state_abb)
    return round(length(findall((biv.==biv_cat) .& (state.==state_abb))) / length(findall(biv.==biv_cat)) * 100, digits=1)
end
pct_cat_state(13,"CA") #e.g.

pct_cat_state(33,"CA") 
pct_cat_state(33,"CO") 
pct_cat_state(33,"MT") 
pct_cat_state(33,"NM") 
pct_cat_state(33,"UT") 
pct_cat_state(33,"WY")
