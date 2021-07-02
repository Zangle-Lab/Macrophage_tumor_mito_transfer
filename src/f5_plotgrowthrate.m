%
% run after tracking script
% bins specific growth rate data into the groups:
% (1) macrophage:GFP 
% (2) control tumor:RFP 
% (3) mito transferred tumor:GRFP

clear all; clc;
FRootO='H:\Data\Minna\MDAMB231MacrophageMitoTransferExp_1April2021\Trial1\AlternateFLTrackingResults\'; % folder tracking results stored
dt=15/60;               % time between phase image frames
tmin = 0;               % time to start binning data (h)
tmax = 48;              % time to stop binning data (h)
minpathlength = 60;     % minimum length to consider track for growth calculation
GWindow = 3;            % half-size of window for gaussian filtering of mass vs time (units: data points)
GSigma = 3;             % standard deviation of window for gaussian filtering of mass vs. time (units: data points)
maxjump = 2;            % max jump in growth rate from frame to frame
RFPInt=0.15;            % minimum RFP intensity to count as tumor cell
GFPInt=0.3;             % minimum GFP intensity to count as macrophage
MitoInt=0.0163;         % min GFP intensity of mito puncta, compare track results to FL images to set
MitoLng=3;              % minimum no of frames mito puncta signal exist to count as transfer, will depend on tagging and how well focussed FL imaging was

%% initiate specific growth rate matrix for 3 groups
sgr_all_GFP=[];         % sp. growth rate of macrophage 
sgr_all_RFP=[];         % sp. growth rate of control tumor cells
sgr_all_GRFP=[];        % sp. growth rate of mito transferred tumor cells

%% load tracking data, calculate sp. growth rate & store
D = load([FRootO, 'data_allframes_RFP']);
[num, indices] = track_numpart(D.tracks,minpathlength); %find tracks >= minpathlength
ss=0;
for ii = 1:num %loop through all tracks meeting minpathlength
    currentnum = D.tracks(indices(ii),5);
    [~,~,z,t,r,g,l] = track_partn(D.tracks,currentnum);       % z:mass, t:time, r:RFP intensity, g:GFP intensity, l:mito puncta intensity
    if mean(nonzeros(r))>=RFPInt && mean(nonzeros(g))<=GFPInt % to avoid macrophages with tumor mito transferred
        if sum(l>=MitoInt)>=MitoLng                           % if mito puncta signal significant for >=MitoLng frames, mito transferred cell
            [fitData,coeff,cnst]= polyfit(t,z,1);             % growth rate is slope of mass vs time plot
            gr=fitData(1)/cnst(2);                            
            sgr=mean(gr)/mean(z);                             % specific growth rate = growth rate / average cell mass
            if sgr>0
                sgr_all_GRFP=[sgr_all_GRFP,sgr];           
            end
        else
            [fitData,coeff,cnst]= polyfit(t,z,1);             % control tumor cell
            gr=fitData(1)/cnst(2);
            sgr=mean(gr)/mean(z);
            if sgr>0
                sgr_all_RFP=[sgr_all_RFP,sgr];
            end
        end
    end
    if mean(nonzeros(g))>=GFPInt && mean(nonzeros(r))<=RFPInt % macrophage growth rate
       [fitData,coeff,cnst]= polyfit(t,z,1);
        gr=fitData(1)/cnst(2);
        sgr=mean(gr)/mean(z);
        sgr_all_GFP=[sgr_all_GFP,sgr]; 
    end
end

%% to plot bar graph on avg growth rate for all three cases
figure(1);
sem_GFP=std(sgr_all_GFP)/sqrt(length(sgr_all_GFP));
sem_RFP=std(sgr_all_RFP)/sqrt(length(sgr_all_RFP));
sem_GRFP=std(sgr_all_GRFP)/sqrt(length(sgr_all_GRFP));
bar([mean(sgr_all_GFP),mean(sgr_all_RFP),mean(sgr_all_GRFP)]);
set(gca, 'XTickLabel', {'Macrophage' 'Tumor Control','Tumor with mito'});
hold on;
errorbar([1,2,3],[mean(sgr_all_GFP),mean(sgr_all_RFP),mean(sgr_all_GRFP)],[sem_GFP,sem_RFP,sem_GRFP],[sem_GFP,sem_RFP,sem_GRFP],'.','Color','b');

%% scatter plot to see the distribution of growth rate values in three groups
figure(2);
DataAlls{1,1}=nonzeros(sgr_all_GFP);
DataAlls{1,2}=nonzeros(sgr_all_RFP);
DataAlls{1,3}=nonzeros(sgr_all_GRFP);

for dd=1:3
    scatter(dd*ones(1,length(DataAlls{1,dd})),DataAlls{1,dd});
    hold on; 
    plot([dd-0.25,dd+0.25], [mean(DataAlls{1,dd}),mean(DataAlls{1,dd})], '-r','LineWidth',2)
    errorbar(dd,mean(DataAlls{1,dd}),std(DataAlls{1,dd})/sqrt(length(DataAlls{1,dd})),std(DataAlls{1,dd})/sqrt(length(DataAlls{1,dd})));

end
xlim([0,4]);
set(gca,'xtick',1:3,'xticklabel',{'Macrophage','Control tumor','Mito tumor'})
    