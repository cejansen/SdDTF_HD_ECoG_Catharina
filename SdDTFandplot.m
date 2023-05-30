eeglab %Opens EEGLAB, filtering and event selection are done in the GUI

%To use this script the SIFT toolbox has to be installed as an EEGlab
%extension
%%  parameters
    EpochTimeLimits     = [0 2];                 % this is the time range (in seconds) to analyze (relative to event at t=0)
    WindowLengthSec     = [2];        % sliding window length in seconds
    WindowStepSizeSec   = [2];      % sliding window step size in seconds
    SelAlg              = {'ARfit'};      %,'Group Lasso (DAL/;;;;;SCSA)','Group Lasso (ADMM)'};        % selection algorithms to determine model order
    NewSamplingRate     = [];                           % new sampling rate (if downsampling)
    GUI_MODE            = 'nogui';                      % whether or not to show the Graphical User Interfaces. Can be 'nogui' or anything else (to show the gui)
    VERBOSITY_LEVEL     = 2;                            % Verbosity Level (0=no/minimal output, 2=graphical output)
    morder=30;
%%
[ALLEEG] = pop_pre_prepData(ALLEEG,GUI_MODE,'VerbosityLevel',VERBOSITY_LEVEL,'SignalType',{'Channels'},'Detrend',[],  ...
        'NormalizeData',[],'resetConfigs',true,'badsegments',[],'newtrials',[],'equalizetrials',false);

%% Estimate MVAR model for the 3 bandwidths
for q=1:3
    switch(q)
        case{1}
            freq=[30:80];
        case{2}
            freq=[80:250];
        case{3}
            freq=[250:500];
    end

ALLEEG(q)= pop_est_fitMVAR(ALLEEG(q),GUI_MODE,'VerbosityLevel',VERBOSITY_LEVEL,'algorithm',SelAlg,'ModelOrder',morder, 'EpochTimeLimits',EpochTimeLimits,'WindowLength',WindowLengthSec,'WindowStepSize',WindowStepSizeSec);%,'EpochTimelimits',EpochTimeRange,'WindowLength',WindowLengthSec,'WindowStepSize',WindowStepSizeSec);
end
%% Calculate connectivity
for j=1:3
    switch(j)
        case{1}
            freq=[30:80];
        case{2}
            freq=[80:250];
        case{3}
            freq=[250:500];
    end
ALLEEG(j) = pop_est_mvarConnectivity(ALLEEG(j),GUI_MODE, ...
                            'connmethods',{'DTF' 'nDTF' 'dDTF' 'dDTF08' 'ffDTF'}, 'absvalsq',true,'spectraldecibels',false,   ...
                            'verb',VERBOSITY_LEVEL,'freqs',freq);
end                      

%% Extract connectivity data and set diagonal to zero
CONNgam=ALLEEG(1).CAT.Conn; 
CONNrip=ALLEEG(2).CAT.Conn;
CONNfr= ALLEEG(3).CAT.Conn; 
CONNgam.ntotDTF = CONNgam.DTF./sum(sum(sum(CONNgam.DTF)));
CONNgam.ntotDTF(1:(size(CONNgam.ntotDTF,1)+1):end) = 0;  
CONNrip.ntotDTF = CONNrip.DTF./sum(sum(sum(CONNrip.DTF)));
CONNrip.ntotDTF(1:(size(CONNrip.ntotDTF,1)+1):end) = 0; 
CONNfr.ntotDTF = CONNfr.DTF./sum(sum(sum(CONNfr.DTF)));
CONNfr.ntotDTF(1:(size(CONNfr.ntotDTF,1)+1):end) = 0; 

%% Calculate in,  out and total flow
fnames = fieldnames(CONNgam);
for m = 1:size(CONNgam,1)
    for k = 5:length(fnames)
        INg(m).(sprintf('%s',fnames{k})) = squeeze(sum(CONNgam(m).(sprintf('%s',fnames{k})),2));
        OUTg(m).(sprintf('%s',fnames{k})) = squeeze(sum(CONNgam(m).(sprintf('%s',fnames{k})),1));
        TOTg(m).(sprintf('%s',fnames{k})) = sum(cat(3,INg(m).(sprintf('%s',fnames{k})),OUTg(m).(sprintf('%s',fnames{k}))),3);
        
        INrip(m).(sprintf('%s',fnames{k})) = squeeze(sum(CONNrip(m).(sprintf('%s',fnames{k})),2));
        OUTrip(m).(sprintf('%s',fnames{k})) = squeeze(sum(CONNrip(m).(sprintf('%s',fnames{k})),1));
        TOTrip(m).(sprintf('%s',fnames{k})) = sum(cat(3,INrip(m).(sprintf('%s',fnames{k})),OUTrip(m).(sprintf('%s',fnames{k}))),3);
        
        INfr(m).(sprintf('%s',fnames{k})) = squeeze(sum(CONNfr(m).(sprintf('%s',fnames{k})),2));
        OUTfr(m).(sprintf('%s',fnames{k})) = squeeze(sum(CONNfr(m).(sprintf('%s',fnames{k})),1));
        TOTfr(m).(sprintf('%s',fnames{k})) = sum(cat(3,INfr(m).(sprintf('%s',fnames{k})),OUTfr(m).(sprintf('%s',fnames{k}))),3);
    end

end
%% average in, out, total flow per bandwidth
for i=1:63
    OUTgs(i)=sum(OUTg.ntotDTF(i,:));
    INgs(i)=sum(INg.ntotDTF(i,:));
    TOTgs(i)=sum(TOTg.ntotDTF(i,:));

    OUTrips(i)=sum(OUTrip.ntotDTF(i,:));
    INgrips(i)=sum(INrip.ntotDTF(i,:));
    TOTrips(i)=sum(TOTrip.ntotDTF(i,:));
    
    OUTfrs(i)=sum(OUTfr.ntotDTF(i,:));
    INfrs(i)=sum(INfr.ntotDTF(i,:));
    TOTfrs(i)=sum(TOTfr.ntotDTF(i,:));
end
save("11a","OUTfrs","TOTfrs","INfrs","OUTrips","TOTrips","INgrips","OUTgs","TOTgs","INgs");
%% Dotplot and load resected electrodes
scale=10; %Scale of the dotplot, 10 in all 4 recordings
%gres=[22 30 37 38 39 45 46 47 53 54];%11b
gres=[12 13 20 21 22 28 29 30 36 37];%11a Resected electrodes
%gres=[12 13 14 15 18 19 20 21 22 23 26 27 28 29 30 35 36 37 38];%21a
%gres=[10 11 12 13 14 15 16 19 20 21 22 23 27 28 29 30 31 36 37 38 44 45];%31a
dotplot(INgs,OUTgs,TOTgs,gres,scale);title("Gamma")
dotplot(INgrips,OUTrips,TOTrips,gres,scale);title("Ripple")
dotplot(INfrs,OUTfrs,TOTfrs,gres,scale);title("Fast ripple")
%% Save resected and non-resected electrodes seperately
ja=1;nr=1;
for r=1:63
    rescheck=any(gres(:) == r);
    if rescheck==1
        res.in.g(ja)=INgs(r);
        res.out.g(ja)=OUTgs(r);
        res.tot.g(ja)=TOTgs(r);

        res.in.rip(ja)=INgrips(r);
        res.out.rip(ja)=OUTrips(r);
        res.tot.rip(ja)=TOTrips(r);

        res.in.fr(ja)=INfrs(r);
        res.out.fr(ja)=OUTfrs(r);
        res.tot.fr(ja)=TOTfrs(r);
        ja=ja+1;
    else
        nres.in.g(nr)=INgs(r);
        nres.out.g(nr)=OUTgs(r);
        nres.tot.g(nr)=TOTgs(r);
    
        nres.in.rip(nr)=INgrips(r);
        nres.out.rip(nr)=OUTrips(r);
        nres.tot.rip(nr)=TOTrips(r);
    
        nres.in.fr(nr)=INfrs(r);
        nres.out.fr(nr)=OUTfrs(r);
        nres.tot.fr(nr)=TOTfrs(r);
        nr=nr+1;

    end
end
save("P1_m1a_res","res")
save("p1_m1a_nres","nres")
%% Pool patients together for boxplot
nonres_in_g=[nonres_in_g nres.in.g];
nonres_in_r=[nonres_in_r nres.in.rip];
nonres_in_fr=[nonres_in_fr nres.in.fr];

nonres_out_g=[nonres_out_g nres.out.g];
nonres_out_r=[nonres_out_r nres.out.rip];
nonres_out_fr=[nonres_out_fr nres.out.fr];

nonres_tot_g=[nonres_tot_g nres.tot.g];
nonres_tot_r=[nonres_tot_r nres.tot.rip];
nonres_tot_fr=[nonres_tot_fr nres.tot.fr];
%%
save("nonres","nonres_tot_fr","nonres_tot_g","nonres_tot_r","nonres_out_fr","nonres_out_r","nonres_out_g","nonres_in_fr","nonres_in_r","nonres_in_g")
%%
res_in_g=[res_in_g res.in.g];
res_in_r=[res_in_r res.in.rip];
res_in_fr=[res_in_fr res.in.fr];

res_out_g=[res_out_g res.out.g];
res_out_r=[res_out_r res.out.rip];
res_out_fr=[res_out_fr res.out.fr];

res_tot_g=[res_tot_g res.tot.g];
res_tot_r=[res_tot_r res.tot.rip];
res_tot_fr=[res_tot_fr res.tot.fr];
%%
save("res","res_tot_fr","res_tot_r","res_tot_g","res_out_fr","res_out_r","res_out_g","res_in_fr","res_in_r","res_in_g")
%% Plot the boxplots
a=191; %Total number of non-resected electrodes
b=61; %Total number of resected electrodes
x=i[nonres_in_g'; res_in_g';nonres_in_r';res_in_r';nonres_in_fr';res_in_fr'];
g1 = repmat({'res. Gamma'},b,1);
g2 = repmat({'non-res. Gamma'},a,1);
g3 = repmat({'res. Ripple'},b,1);
g4 = repmat({'non-res. Ripple'},a,1);
g5 = repmat({'res. Fast Ripple'},b,1);
g6 = repmat({'non-res. Fast Ripple'},a,1);
g = [g1;g2; g3;g4;g5;g6];
boxplot(xi,g)
title('Inflow')
figure(5)
xo=[nonres_out_g'; res_out_g';nonres_out_r';res_out_r';nonres_out_fr';res_out_fr'];
boxplot(xo,g)
title('Outflow')
figure(6)
xt=[nonres_tot_g'; res_tot_g';nonres_tot_r';res_tot_r';nonres_tot_fr';res_tot_fr'];
boxplot(xt,g)
title('Total flow')