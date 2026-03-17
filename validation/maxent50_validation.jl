# Mitchell Chandler
# mitchell_chandler@tws.org
# Last updated: 06-March-2026
##################################################

#import Pkg

using NetCDF
using CairoMakie
using NumericalIntegration
using Statistics
using CSV
using DataFrames
using HypothesisTests 

cd("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/_Data/scripts/50pct_witheld")

#CA; CO; MT; NM; UT; WY

prob_threshold = [1:-0.05:0;];

##
# SPATIAL AUTOCORRELATION OF RESIDUALS #

#Distance [m]
dist = CSV.read("CA50_sa.csv",DataFrame)[:,:Distance];  #distances should be the same for all states

#Moran's I
I_CA = CSV.read("CA50_sa.csv",DataFrame)[:,:MoransI]; #CA  
I_CO = CSV.read("CO50_sa.csv",DataFrame)[:,:MoransI]; #CO  
I_MT = CSV.read("MT50_sa.csv",DataFrame)[:,:MoransI]; #MT  
I_NM = CSV.read("NM50_sa.csv",DataFrame)[:,:MoransI]; #NM  
I_UT = CSV.read("UT50_sa.csv",DataFrame)[:,:MoransI]; #UT  
I_WY = CSV.read("WY50_sa.csv",DataFrame)[:,:MoransI]; #WY  

# Plot
f = Figure(fontsize=16);
ax = Axis(f[1, 1],xlabel="Distance [km]",ylabel="Moran's I", 
    xminorticks=IntervalsBetween(5),xminorticksvisible=false,xgridvisible=false,ygridvisible=false,xminorgridvisible=false)
hlines!(0,color=:black,linewidth=1)
lines!(dist/1000,I_CA,linewidth=4,label="California",color=Makie.to_colormap(:Dark2_7)[1])
lines!(dist/1000,I_CO,linewidth=4,label="Colorado",color=Makie.to_colormap(:Dark2_7)[2])
lines!(dist/1000,I_MT,linewidth=4,label="Montana",color=Makie.to_colormap(:Dark2_7)[3])
lines!(dist/1000,I_NM,linewidth=4,label="New Mexico",color=Makie.to_colormap(:Dark2_7)[5])
lines!(dist/1000,I_UT,linewidth=4,label="Utah",color=Makie.to_colormap(:Dark2_7)[6])
lines!(dist/1000,I_WY,linewidth=4,label="Wyoming",color=Makie.to_colormap(:Dark2_7)[7])
#axis:
xlims!(0,300)
ylims!(-0.1,1)
ax.yticks=-0.1:0.1:1;
ax.xticks=0:50:300;
#legend:
f[1,2] = Legend(f,ax,framevisible=true,backgroundcolor=:white,framecolor=:white)
display(f)
#save("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/Figures/spatial_autocorrelogram-2026.03.06.pdf",f)

#find minimum distance that spatial autocorrelation for all states is less than or equal to a given Moran's I
I = 0.1;
minimum(dist[max.(I_CA,I_CO,I_MT,I_NM,I_UT,I_WY).<=I])/1000


##
# OMISSION RATE AND AUC OF TRAINING AND TESTING DATA #

function otb(prob,omission,bground_prob,threshold)
    #one-tailed binomial test to assess significance of omission rate - Phillips et al. 2006, sec 2.6.1
    t = length(prob); #test localities
    r = omission; #omission rate
    a = length(findall([prob; bground_prob].>=threshold)) / (length(prob) + length(bground_prob)); #proportional predicted area
    p = pvalue(BinomialTest(round(t*(1-r)),t,a),tail=:right); #right-tailed bc the alt hypothesis is that success rate is greater than
    return p
end

function model_val(train_nc,test_nc,bground_nc,prob_threshold,title_state)
    # DATA #
    #OGDP at grid well training locations
    train_prob = ncread(train_nc,"grid_code");
    trainx = ncread(train_nc,"x");
    trainy = ncread(train_nc,"y");
    #OGDP at grid well testing locations
    test_prob = ncread(test_nc,"PROB");
    testx = ncread(test_nc,"x");
    testy = ncread(test_nc,"y");
    #OGDP at background locations
    bground_prob = ncread(bground_nc,"grid_code");
    bgroundx = ncread(bground_nc,"x");
    bgroundy = ncread(bground_nc,"y");
    #plot
    f = Figure();
    ax = Axis(f[1, 1],aspect = DataAspect())
    sc = scatter!(bgroundx,bgroundy,color=bground_prob,colorrange=(0,1),colormap=:BuPu,markersize=6)
    sc = scatter!(testx,testy,color=test_prob,colorrange=(0,1),colormap=:BuPu,markersize=6)
    sc = scatter!(trainx,trainy,color=train_prob,colorrange=(0,1),colormap=:BuPu,markersize=6)
    Colorbar(f[1,2],sc,ticks=([0:0.1:1;]),height=300,width=20,tellheight=false)
    hidedecorations!(ax)
    display(f)

    # VALIDATION STATISTICS #
    #initialise
    train_omission = NaN*prob_threshold;
    test_omission = NaN*prob_threshold;
    bground_rate = NaN*prob_threshold;
    train_p = NaN*prob_threshold;
    test_p = NaN*prob_threshold;
    #compute omission rate and background occurrence rate
    for i in eachindex(prob_threshold)
        train_omission[i] = length(findall(train_prob.<prob_threshold[i]))/length(train_prob);
        
        train_p[i] = otb(train_prob,train_omission[i],bground_prob,prob_threshold[i]);

        test_omission[i] = length(findall(test_prob.<prob_threshold[i]))/length(test_prob);
        
        test_p[i] = otb(test_prob,test_omission[i],bground_prob,prob_threshold[i]);

        bground_rate[i] = length(findall(bground_prob.>=prob_threshold[i]))/length(bground_prob);    
    end
    #plot 
    f = Figure(fontsize=16);
    
    #Omission Rate
    ax1 = Axis(f[1, 1],xlabel="Probability Threshold",ylabel="Omission Rate",title=title_state)
    #training data:
    lines!(prob_threshold,train_omission,linewidth=4,color=Makie.to_colormap(:batlowS)[1],label="Training")
    scatterlines!(prob_threshold[train_p.<=0.05],train_omission[train_p.<=0.05],linewidth=0,color=Makie.to_colormap(:batlowS)[1],markersize=12) #significant (p<0.05)
    scatterlines!(prob_threshold[train_p.>0.05],train_omission[train_p.>0.05],linewidth=0,color=:grey60,markersize=12) #non-significant (p>0.05)
    #testing data:
    lines!(prob_threshold,test_omission,linewidth=4,color=Makie.to_colormap(:batlowS)[2],label="Testing")
    scatterlines!(prob_threshold[test_p.<=0.05],test_omission[test_p.<=0.05],linewidth=0,color=Makie.to_colormap(:batlowS)[2],markersize=12) #significant (p<0.05)
    scatterlines!(prob_threshold[test_p.>0.05],test_omission[test_p.>0.05],linewidth=0,color=:grey60,markersize=12) #non-significant (p>0.05)
    #axis:
    xlims!(0,1)
    ylims!(0,1)
    ax1.xticks=0:0.25:1;
    ax1.yticks=0:0.1:1;
    #legend:
    axislegend(framevisible=true,position=:lt,backgroundcolor=:white,framecolor=:white)

    #ROC Curve
    ax2 = Axis(f[1, 2],xlabel="Fraction of background\n points classified as presence\n (1 - Specificity)",ylabel="True positive rate (Sensitivity)",title=title_state)
    lines!(bground_rate,(1 .- train_omission),linewidth=4,color=Makie.to_colormap(:batlowS)[1],label="Training")
    lines!(bground_rate,(1 .- test_omission),linewidth=4,color=Makie.to_colormap(:batlowS)[2],label="Testing")
    #axis:
    xlims!(0,1)
    ylims!(0,1)
    ax2.xticks=0:0.2:1;
    ax2.yticks=0:0.1:1;
    #legend:
    axislegend(framevisible=true,position=:rb,backgroundcolor=:white,framecolor=:white)
    
    display(f)

    #compute AUC
    train_AUC = integrate(bground_rate,(1 .- train_omission))
    test_AUC = integrate(bground_rate,(1 .- test_omission))

    return train_AUC, test_AUC, test_omission, test_p
end

#CA
CA_train_AUC, CA_test_AUC, CA_test_omission, CA_test_p = model_val("CA50_training_points.nc","CA50_prediction_points.nc","CA50_background_points.nc",prob_threshold,"CA");
CA_train_AUC, CA_test_AUC

#CO
CO_train_AUC, CO_test_AUC, CO_test_omission, CO_test_p = model_val("CO50_training_points.nc","CO50_prediction_points.nc","CO50_background_points.nc",prob_threshold,"CO");
CO_train_AUC, CO_test_AUC

#MT
MT_train_AUC, MT_test_AUC, MT_test_omission, MT_test_p = model_val("MT50_training_points.nc","MT50_prediction_points.nc","MT50_background_points.nc",prob_threshold,"MT");
MT_train_AUC, MT_test_AUC

#NM
NM_train_AUC, NM_test_AUC, NM_test_omission, NM_test_p = model_val("NM50_training_points.nc","NM50_prediction_points.nc","NM50_background_points.nc",prob_threshold,"NM");
NM_train_AUC, NM_test_AUC

#UT
UT_train_AUC, UT_test_AUC, UT_test_omission, UT_test_p = model_val("UT50_training_points.nc","UT50_prediction_points.nc","UT50_background_points.nc",prob_threshold,"UT");
UT_train_AUC, UT_test_AUC

#WY
WY_train_AUC, WY_test_AUC, WY_test_omission, WY_test_p = model_val("WY50_training_points.nc","WY50_prediction_points.nc","WY50_background_points.nc",prob_threshold,"WY");
WY_train_AUC, WY_test_AUC

#round(WY_test_omission[prob_threshold.==0.5][1],digits=2)

# Output AUC for training and testing data
AUC_df = DataFrame(
    State=["CA","CO","MT","NM","UT","WY"],
    Train_AUC=round.([CA_train_AUC,CO_train_AUC,MT_train_AUC,NM_train_AUC,UT_train_AUC,WY_train_AUC],digits=4),
    Test_AUC=round.([CA_test_AUC,CO_test_AUC,MT_test_AUC,NM_test_AUC,UT_test_AUC,WY_test_AUC],digits=4),
    AUC_diff=round.([CA_train_AUC,CO_train_AUC,MT_train_AUC,NM_train_AUC,UT_train_AUC,WY_train_AUC] - [CA_test_AUC,CO_test_AUC,MT_test_AUC,NM_test_AUC,UT_test_AUC,WY_test_AUC],digits=4)
)
#clipboard(sprint(show,"text/tab-separated-values",AUC_df))

##
# PLOT OMISSION RATES OF TESTING DATA #

function plot_state_omission(prob,state,omission,p,cmap,txt)
    lines!(prob,omission,linewidth=4,label=state,color=cmap)
    #scatterlines!(prob[p.<=0.05],omission[p.<=0.05],linewidth=0,color=cmap,markersize=12)
    scatterlines!(prob[p.>0.05],omission[p.>0.05],linewidth=0,color=:black,markersize=9,strokecolor=cmap,strokewidth=1)
    #text:
    text!(0.05,0.9,text=txt,color=cmap,align=(:left,:center))
    #axis:
    xlims!(0,1)
    ylims!(0,1)
end

f = Figure();  
#CA
ax1 = Axis(f[1, 1],ylabel="Omission Rate",
    xminorticks=IntervalsBetween(5),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true,
    xgridvisible=true,ygridvisible=true,
    xminorgridvisible=true,yminorgridvisible=true,
    yticksmirrored=true)
plot_state_omission(prob_threshold,"California",CA_test_omission,CA_test_p,Makie.to_colormap(:Dark2_7)[1],"(a) CA")
#axis:
ax1.xticks=0:0.5:1;
ax1.yticks=0:0.2:1;

#CO
ax2 = Axis(f[1, 2],
    xminorticks=IntervalsBetween(5),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true,
    xgridvisible=true,ygridvisible=true,
    xminorgridvisible=true,yminorgridvisible=true,
    yticklabelsvisible=false,
    yticksmirrored=true)
plot_state_omission(prob_threshold,"Colorado",CO_test_omission,CO_test_p,Makie.to_colormap(:Dark2_7)[2],"(b) CO")
#axis:
ax2.xticks=0:0.5:1;
ax2.yticks=0:0.2:1;

#MT
ax3 = Axis(f[1, 3],
    xminorticks=IntervalsBetween(5),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true,
    xgridvisible=true,ygridvisible=true,
    xminorgridvisible=true,yminorgridvisible=true,
    yticklabelsvisible=false,
    yticksmirrored=false)
plot_state_omission(prob_threshold,"Montana",MT_test_omission,MT_test_p,Makie.to_colormap(:Dark2_7)[3],"(c) MT")
#axis:
ax3.xticks=0:0.5:1;
ax3.yticks=0:0.2:1;

#NM
ax5 = Axis(f[2, 1],ylabel="Omission Rate",
    xlabel="Probability",
    xminorticks=IntervalsBetween(5),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true,
    xgridvisible=true,ygridvisible=true,
    xminorgridvisible=true,yminorgridvisible=true,
    yticksmirrored=true)
plot_state_omission(prob_threshold,"New Mexico",NM_test_omission,NM_test_p,Makie.to_colormap(:Dark2_7)[5],"(d) NM")
#axis:
ax5.xticks=0:0.5:1;
ax5.yticks=0:0.2:1;

#UT
ax6 = Axis(f[2, 2],xlabel="Probability",
    xminorticks=IntervalsBetween(5),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true,
    xgridvisible=true,ygridvisible=true,
    xminorgridvisible=true,yminorgridvisible=true,
    yticklabelsvisible=false,
    yticksmirrored=true)
plot_state_omission(prob_threshold,"Utah",UT_test_omission,UT_test_p,Makie.to_colormap(:Dark2_7)[6],"(e) UT")
#axis:
ax6.xticks=0:0.5:1;
ax6.yticks=0:0.2:1;

#WY
ax7 = Axis(f[2, 3],xlabel="Probability",
    xminorticks=IntervalsBetween(5),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true,
    xgridvisible=true,ygridvisible=true,
    xminorgridvisible=true,yminorgridvisible=true,
    yticklabelsvisible=false,
    yticksmirrored=false)
plot_state_omission(prob_threshold,"Wyoming",WY_test_omission,WY_test_p,Makie.to_colormap(:Dark2_7)[7],"(f) WY")
#axis:
ax7.xticks=0:0.5:1;
ax7.yticks=0:0.2:1;

display(f)
#save("C:/Users/MitchellChandler/OneDrive - THE WILDERNESS SOCIETY/OGDP model/Figures/testing_omission-2026.01.05.pdf",f)

##
# OMISSION RATE AND AUC OF OGDP COMPOSITE #

function composite_validation(well_nc,background_nc,prob_threshold)
    # DATA #
    #OGDP at well locations
    pred_prob = ncread(well_nc,"grid_code");
    predx = ncread(well_nc,"x");
    predy = ncread(well_nc,"y");
    #OGDP at background locations
    bground_prob = ncread(background_nc,"grid_code");
    bgroundx = ncread(background_nc,"x");
    bgroundy = ncread(background_nc,"y");
    #plot
    f = Figure();
    ax = Axis(f[1, 1],aspect = DataAspect())
    sc = scatter!(bgroundx,bgroundy,color=bground_prob,colorrange=(0,1),colormap=:YlOrBr,markersize=6)
    sc = scatter!(predx,predy,color=pred_prob,colorrange=(0,1),colormap=:YlOrBr,markersize=6)
    Colorbar(f[1,2],sc,ticks=([0:0.1:1;]),height=300,width=20,tellheight=false)
    hidedecorations!(ax)
    display(f)

    # VALIDATION STATISTICS #
    #initialise
    omission_rate = NaN*prob_threshold;
    background_rate = NaN*prob_threshold;
    p = NaN*prob_threshold;
    #compute omission rate and background occurrence rate
    for i in eachindex(prob_threshold)
        omission_rate[i] = length(findall(pred_prob.<prob_threshold[i]))/length(pred_prob);

        p[i] = otb(pred_prob,omission_rate[i],bground_prob,prob_threshold[i]);

        background_rate[i] = length(findall(bground_prob.>=prob_threshold[i]))/length(bground_prob);
    end
    #plot 
    f = Figure(fontsize=16);
    #Omission Rate
    ax1 = Axis(f[1, 1],xlabel="Probability Threshold",ylabel="Omission Rate")
    lines!(prob_threshold,omission_rate,linewidth=4,color=:black)
    scatterlines!(prob_threshold[p.<=0.05],omission_rate[p.<=0.05],linewidth=0,color=:black,markersize=12) #significant (p<0.05)
    scatterlines!(prob_threshold[p.>0.05],omission_rate[p.>0.05],linewidth=0,color=:grey60,markersize=12) #non-significant (p>0.05)
    #axis:
    xlims!(0,1)
    ylims!(0,1)
    ax1.xticks=0:0.25:1;
    ax1.yticks=0:0.1:1;
    #ROC Curve
    ax2 = Axis(f[1, 2],xlabel="Fraction of background\n points classified as presence\n (1 - Specificity)",ylabel="True positive rate (Sensitivity)")
    lines!(background_rate,(1 .- omission_rate),linewidth=4,color=:black)
    #axis:
    xlims!(0,1)
    ylims!(0,1)
    ax2.xticks=0:0.2:1;
    ax2.yticks=0:0.1:1;
    display(f)

    #compute AUC
    AUC = integrate(background_rate,(1 .- omission_rate));

    return AUC, omission_rate, p, background_rate
end

composite_AUC, composite_omission_rate, composite_p = composite_validation("ogdp50_composite_prediction_points.nc","OGDP50_composite_background_points.nc",prob_threshold);
composite_AUC

#=
# Plot omission rate only
f = Figure();  
ax1 = Axis(f[1, 2],ylabel="Omission Rate",xlabel="Probability",
    xminorticks=IntervalsBetween(2),xminorticksvisible=true,
    yminorticks=IntervalsBetween(2),yminorticksvisible=true,
    xgridvisible=true,ygridvisible=true,
    xminorgridvisible=true,yminorgridvisible=true)
plot_state_omission(prob_threshold,"",composite_omission_rate,composite_p,:Grey50,"")
#axis:
ax1.xticks=0:0.2:1;
ax1.yticks=0:0.2:1;
display(f)
=#
