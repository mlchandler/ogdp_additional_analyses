# Mitchell Chandler
# mitchell_chandler@tws.org
# Last updated: 06-March-2026
##################################################

using CairoMakie
using CSV
using DataFrames

cd("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/_Data/scripts/50pct_witheld")

##################################################
## READ DATA ##

function read_rc(fname)
    df = CSV.read(fname,DataFrame);

    mag_val = disallowmissing(df[:,:VALUE][findall(df[:,:VARIABLE] .== "MAGNETIC_ALBERS_10K")]);
    mag_prob = disallowmissing(df[:,:PROB][findall(df[:,:VARIABLE] .== "MAGNETIC_ALBERS_10K")]);

    iso_val = disallowmissing(df[:,:VALUE][findall(df[:,:VARIABLE] .== "ISOSTATIC_ALBERS_10K")]);
    iso_prob = disallowmissing(df[:,:PROB][findall(df[:,:VARIABLE] .== "ISOSTATIC_ALBERS_10K")]);

    bou_val = disallowmissing(df[:,:VALUE][findall(df[:,:VARIABLE] .== "BOUGUER_ALBERS_10K")]);
    bou_prob = disallowmissing(df[:,:PROB][findall(df[:,:VARIABLE] .== "BOUGUER_ALBERS_10K")]);

    bed_val = disallowmissing(df[:,:VALUE][findall(df[:,:VARIABLE] .== "BEDROCK_ALBERS_10K")]);
    bed_prob = disallowmissing(df[:,:PROB][findall(df[:,:VARIABLE] .== "BEDROCK_ALBERS_10K")]);
    #bedrock depth can't be less than 0
    bed_idx = findall(bed_val .>= 0);
    bed_val = bed_val[bed_idx];
    bed_prob = bed_prob[bed_idx];

    top_val = disallowmissing(df[:,:VALUE][findall(df[:,:VARIABLE] .== "ETOPO2022_ALBERS_10K")]);
    top_prob = disallowmissing(df[:,:PROB][findall(df[:,:VARIABLE] .== "ETOPO2022_ALBERS_10K")]);

    lit_val = disallowmissing(df[:,:CATEGORY][findall(df[:,:VARIABLE] .== "LITH_MASKED")]);
    lit_prob = disallowmissing(df[:,:PROB][findall(df[:,:VARIABLE] .== "LITH_MASKED")]);

    return hcat(mag_val,mag_prob), hcat(iso_val,iso_prob), hcat(bou_val,bou_prob), hcat(bed_val,bed_prob), hcat(top_val,top_prob), hcat(lit_val,lit_prob)
end

CA_mag, CA_iso, CA_bou, CA_bed, CA_top, CA_lit = read_rc("CA50_rc.csv"); #CA
CO_mag, CO_iso, CO_bou, CO_bed, CO_top, CO_lit = read_rc("CO50_rc.csv"); #CO
MT_mag, MT_iso, MT_bou, MT_bed, MT_top, MT_lit = read_rc("MT50_rc.csv"); #MT
NM_mag, NM_iso, NM_bou, NM_bed, NM_top, NM_lit = read_rc("NM50_rc.csv"); #NM
UT_mag, UT_iso, UT_bou, UT_bed, UT_top, UT_lit = read_rc("UT50_rc.csv"); #UT
WY_mag, WY_iso, WY_bou, WY_bed, WY_top, WY_lit = read_rc("WY50_rc.csv"); #WY

##################################################
## PLOT CONTINUOUS ##

function plot_rc(CA,CO,MT,NM,UT,WY)
    #response curves:
    lines!(CA[:,1],CA[:,2],linewidth=3,label="California",color=Makie.to_colormap(:Dark2_7)[1])
    lines!(CO[:,1],CO[:,2],linewidth=3,label="Colorado",color=Makie.to_colormap(:Dark2_7)[2])
    lines!(MT[:,1],MT[:,2],linewidth=3,label="Montana",color=Makie.to_colormap(:Dark2_7)[3])
    lines!(NM[:,1],NM[:,2],linewidth=3,label="New Mexico",color=Makie.to_colormap(:Dark2_7)[5])
    lines!(UT[:,1],UT[:,2],linewidth=3,label="Utah",color=Makie.to_colormap(:Dark2_7)[6])
    lines!(WY[:,1],WY[:,2],linewidth=3,label="Wyoming",color=Makie.to_colormap(:Dark2_7)[7])
    #axis:
    ylims!(0,1)
    ax.yticks=0:0.2:1;
end

fsize=11;
tsize=11;

f = Figure(fontsize=fsize);
    #Magnetic anomalies
ax = Axis(f[1, 1],ylabel="Probability of Presence",xlabel="[nT]",
    xgridvisible=false,ygridvisible=false,
    yticksmirrored=true,
    title="(a) Magnetic anomalies",titlefont=:regular,titlealign=:right,titlegap=0,titlesize=tsize)
plot_rc(CA_mag,CO_mag,MT_mag,NM_mag,UT_mag,WY_mag)
    #Isostatic
ax = Axis(f[1, 2],xlabel=rich("[10",superscript("-5")," m/s",superscript("2"),"]"),
    xgridvisible=false,ygridvisible=false,
    yticksmirrored=true,yticklabelsvisible=false,
    title="(b) Isostatic gravity anomalies",titlefont=:regular,titlealign=:right,titlegap=0,titlesize=tsize)
plot_rc(CA_iso,CO_iso,MT_iso,NM_iso,UT_iso,WY_iso)
    #Bouguer
ax = Axis(f[1, 3],xlabel=rich("[10",superscript("-5")," m/s",superscript("2"),"]"),
    xgridvisible=false,ygridvisible=false,
    yticklabelsvisible=false,
    title="(c) Bouguer gravity anomalies",titlefont=:regular,titlealign=:right,titlegap=0,titlesize=tsize)
plot_rc(CA_bou,CO_bou,MT_bou,NM_bou,UT_bou,WY_bou)
    #Bedrock
ax = Axis(f[2, 1],ylabel="Probability of Presence",xlabel="[cm]",
    xgridvisible=false,ygridvisible=false,
    yticksmirrored=true,
    xscale=log10,xminorticks=IntervalsBetween(10),xminorticksvisible=true,
    title="(d) Bedrock depth",titlefont=:regular,titlealign=:right,titlegap=0,titlesize=tsize)
plot_rc(CA_bed,CO_bed,MT_bed,NM_bed,UT_bed,WY_bed)
    #Topography
ax = Axis(f[2, 2],xlabel="[m]",
    xgridvisible=false,ygridvisible=false,
    yticklabelsvisible=false,
    title="(e) Elevation",titlefont=:regular,titlealign=:right,titlegap=0,titlesize=tsize)
plot_rc(CA_top,CO_top,MT_top,NM_top,UT_top,WY_top)
#legend:
f[2,3] = Legend(f,ax,framevisible=true,backgroundcolor=:white,framecolor=:white,tellwidth=false)
display(f)
#save("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/Figures/continuous_response_curve-2026.03.06.pdf",f)

##################################################
## PLOT CATEGORICAL ##

num_state=6;

#Combine lithology response curves for all States into a single matrix by matching lithology categories
all_lit = vcat(CA_lit,CO_lit,MT_lit,NM_lit,UT_lit,WY_lit);
lit_names = sort(unique(all_lit[:,1]));
lit_prob = zeros(length(lit_names),num_state);

function add_lit(state_lit,state_idx)
#use the following for state_idx: CA (1), CO (2), MT (3), NM (4), UT (5), WY (6)
    for i in eachindex(state_lit[:,1])
        idx = findall(state_lit[i,1] .== lit_names)[1];
        lit_prob[idx,state_idx] = state_lit[i,2];
    end
end

add_lit(CA_lit,1)
add_lit(CO_lit,2)
add_lit(MT_lit,3)
add_lit(NM_lit,4)
add_lit(UT_lit,5)
add_lit(WY_lit,6)

#Break out matrix into series of vectors for plotting
lit = (cat = repeat(1:length(lit_names),outer=num_state),
    state = repeat(1:num_state,inner=length(lit_names)),
    prob = lit_prob[:],
    colour = repeat([1,2,3,5,6,7],inner=length(lit_names))
);

#Plot
f = Figure(fontsize=fsize);
ax = Axis(f[1, 1],ylabel="Probability of Presence",
    xticks=(1:length(lit_names),lit_names),
    xticklabelrotation=-pi/5,
    xgridvisible=false,ygridvisible=false)
barplot!(lit.cat,lit.prob,dodge=lit.state,color=Makie.to_colormap(:Dark2_7)[lit.colour])
vlines!(collect(0.5:length(lit_names)+1),linewidth=1,linestyle=:dot,color=:grey50)
    #axis:
ylims!(0,1)
ax.yticks=0:0.1:1;
xlims!(0.5,length(lit_names)+0.5)
    #legend:
labels = ["California","Colorado","Montana","New Mexico","Utah","Wyoming"];
elements = [PolyElement(polycolor = Makie.to_colormap(:Dark2_7)[i]) for i in [1,2,3,5,6,7]]
Legend(f[1,2],elements,labels,framevisible=false)
display(f)
#save("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/Figures/categorical_response_curve-2026.01.05.pdf",f)
