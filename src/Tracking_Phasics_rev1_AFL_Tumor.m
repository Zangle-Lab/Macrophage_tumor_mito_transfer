
%script to run cell tracking code
%TAZ 5/18/10, modified 7/21/10 to automatically detect file names and to
%remove drift due to the entire frame moving
%first section: define image processing parameters (file locations,
%processing options, etc.)
%second section: load images and detect cells
%third section: track cells
%SP modified on 5th Aug 2019 to accomodate 4-D tracking using fluorescent (FL) marker data
%Code reads QPI frames at 15 minute interval & red and green FL images at alternate 15 min intervals
%FL intensity used to avoid tracking untagged cells
%SP modified on 12/10/2020 to bin cells based on FL after tracking

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set up file names of images to be analyzed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
froot='H:\Data\Minna\MDAMB231MacrophageMitoTransferExp_1April2021\Trial1\'; % folder with .mat phase images from Phasics software
fstart='QPM40X'; % phase image name prefix
savefolder = 'H:\Data\Minna\MDAMB231MacrophageMitoTransferExp_1April2021\Trial1\AlternateFLTrackingResults\'; % folder to save track data

%%% First, define which files the script will work on
fext = '.mat'; %file extension
overwrite = 1; %set to 1 to enable overwrite of pre-stored data files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second, define image analysis parameters

%%% define min and max area and mean intensity (MI) thresholds
%%% only objects which fall between these values will be counted as "cells"
%%% These parameters should be adjusted for each sample to capture the
%%% objects of interest. See Figure 11 to evaluate where these values fall
%%% relative to the properties of the image. See figures 12 and 13 to see
%%% which objects in the first and last frames are counted as "cells"
minAreathresh = 600;
maxAreathresh = 8000;
minMIthresh = 50;
maxMIthresh = 1000;

%%% define wavelength of illumination source and pixelsize of phase image
%%% Wavelength and pixel size are used for dry mass measurement from phase image
%%% Pixel size from phase image of caibrated ruler by same Phasics camera
wavelength = 623; % nm
pxlsize = 0.6e-3; % mm/pixel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Third, define tracking software parameters

minpathlength = 10; %min path length to use in plotting results. only paths
%                    of this length or longer will be displayed. this does
%                    not affect the tracking software (tracks shorter than
%                    minpathlength will still be computed and stored)

massfact = 0.8;   %factor to multiply mass by in tracking step. use this to
%account for differences in how much the cell moves vs. how
%much mass changes over time (was set to 2.5 previously, 2/25/15)

%%% tracking parameters below affect the tracking software itself
max_disp = 80; %max displacement for particle tracking
%               max_disp is an estimate of the maximum distance that a
%               particle would move in a single time interval. It should be
%               set to a value somewhat less than the mean spacing between
%               the particles
param.mem = 3; %this is the number of time steps that a particle can be
%               'lost' and then recovered again.  If the particle reappears
%               after this number of frames has elapsed, it will be
%               tracked as a new particle. The default setting is zero.
%               this is useful if particles occasionally 'drop out' of
%               the data.
param.dim = 3; %number of dimensions of coordinate data to track. If you
%               set this value to 2, then it will track based on position
%               in x and y. If param.dim = 3, then the software will track
%               based on mass as well.
param.good = 0; %set this keyword to eliminate all trajectories with
%                fewer than param.good valid positions.  This is useful
%                for eliminating very short, mostly 'lost' trajectories
%                due to blinking 'noise' particles in the data stream.
param.quiet = 1; %set this keyword to 1 if you don't want any text
%                 displayed while the tracking algorithm is running

%% section 2: Initiate tracking by analysis of first image at time=0 min
% Visually comparing background corrected phase image D1 & its label L verifies segmentation
[LocList, numLoc] = getloclist(froot, fstart, fext); % number of positions imaged (numLoc)

%pre-processing and variable initialization before loop begins:
Loc = 1;
filelist = dir([froot, 'QPM40X_', char(LocList(Loc)), '_*', fext]);
fileNames = char(sort_nat({filelist.name}'));
numf = length(fileNames);
fname = strtrim([froot, fileNames(1,:)]);
fname0=sprintf('RFP40X_%d_frame_%d.tif',Loc,1);
RFP=imread([froot fname0]);

% based on a low-high level threshold to enhance the FL image:
limR=single([min(min(RFP)),min(min(RFP))+300]);  % Varies with imaging sets & tagging variability
[D1,L,taglist] = LoadSegment_RFP_AlFL(fname, wavelength,RFP,limR); %D1 is background corrected image, 
% L is label matrix & taglist logs RFP intensity of each cell from FL image
time0 = LoadTime(fname); % image start time point, used as time=0 reference for tracking
D1s = zeros([size(D1),numLoc], 'single'); % zero matrix to store background corrected phase images with tracking

%% section 3: segment and track cells

%loop over all locations
for Loc = 1:numLoc %can parfor
    if ~exist([savefolder, 'Tdata_RFP', num2str(Loc), '.mat']) || overwrite

        filelist = dir([froot, 'QPM40X_', char(LocList(Loc)), '_*', fext]);
        fileNames = char(sort_nat({filelist.name}'));
        numf = length(fileNames(:,1)); % number of frames at each location imaged
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% grab first frame and corresponding FL image for analysis & detection of the correct cell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fname = strtrim([froot, fileNames(1,:)]);
        fname0=sprintf('RFP40X_%d_frame_%d.tif',Loc,1);
        RFP=imread([froot fname0]);
        limR=single([min(min(RFP)),min(min(RFP))+300]); % based on a low level threshold to enhance the image
        [D1,L,taglist] = LoadSegment_RFP_AlFL(fname, wavelength,RFP,limR);
        
        %preallocate variables for speed
        yshift_store = zeros(numf); 
        xshift_store = zeros(numf);
        t_stored = zeros(numf); % time 
        
        D_stored = zeros([size(D1),numf], 'single'); % background corrected phase image
        L_stored = zeros([size(D1),numf], 'uint16'); % label matrix for each frame
        
        
        D1s(:,:,Loc) = single(D1);
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% loop through first numf file names stored in fnum and store analysis
        %%% results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tt = 1; %initialize tt, the index of the tracking array
        T_array = [];
        
        xshift_old = 0;
        yshift_old = 0;
        D_old = D1;
        
        for jj = 1:numf
            
            fname = strtrim([froot fileNames(jj,:)]);
            disp(strtrim(fileNames(jj,:)))
            % lists RFP & GFP FL intensity from alternate images corresponding to phase image
            % Transferred mito intensity logged separately in Gtaglist
            if rem(jj,2)==1
                fname0=sprintf('RFP40X_%d_frame_%d.tif',Loc,jj);
                RFP=imread([froot fname0]);
                limR=single([min(min(RFP)),min(min(RFP))+300]); % based on a low level threshold to enhance the image
                [D,L,Rtaglist] = LoadSegment_RFP_AlFL(fname, wavelength,RFP,limR); % Rtaglist is list of RFP intensity of tumor cells
                Gtaglist=zeros(length(Rtaglist(:,1)),2); Gtaglist(:,1)=Rtaglist(:,1);
                listgfp=zeros(length(Rtaglist(:,1)),2); listgfp(:,1)=Gtaglist(:,1); % GFP intensity data unavialable at the time point
            else
                fname1=sprintf('GFP40X_%d_frame_%d.tif',Loc,jj);
                GFP=imread([froot fname1]);                     % read GFP FL signal in macrophages & mitos transferred to tumor cells
                limG=single([min(min(GFP)),min(min(GFP))+300]);
                [D, L,Gtaglist,listgfp] = LoadSegment_GFP_AlFL(fname, wavelength,GFP, limG); % Macrophage emerald tag intensity (Gtaglist)
                Rtaglist=zeros(length(Gtaglist(:,1)),2); Rtaglist(:,1)=Gtaglist(:,1); % RFP intensity unavailable
            end
            
            timen = LoadTime(fname);
            time = (datenum(timen)-datenum(time0)).*24; %store time in hours
            
            [V, M, A, MI, P, SF] = imageprops_SF(L, D, pxlsize); %compute image properties based on the regions stored in L
            if ~isempty(V)
                listlabel=zeros(1,length(V));
                Rlistlab=zeros(1,length(V));
                Glistlab=zeros(1,length(V));
                if sum(sum(listgfp))~=0
                    for yy=1:length(listgfp(:,2))
                        listlabel(listgfp(yy,1))=listgfp(yy,2);
                        Rlistlab(Rtaglist(yy,1))=Rtaglist(yy,2);
                        Glistlab(Gtaglist(yy,1))=Gtaglist(yy,2);
                    end
                end
            end
            if std(D(:))~=0 %skip if blank image
                [yshift, xshift] = CorrShift(D_old,D); %find average shift between current frame and first frame
                yshift = yshift+yshift_old;
                xshift = xshift+xshift_old;
                D_old = D;
                xshift_old = xshift;
                yshift_old = yshift;
                
                yshift_store(jj,Loc) = yshift;
                xshift_store(jj,Loc) = xshift;
                D_stored(:,:,jj) = single(D(:,:));
                L_stored(:,:,jj) = uint16(L(:,:));
                t_stored(jj,Loc) = time;
                
                %next, loop through all items identified in V and find only the ones
                %which meet area and mean intensity requirements
                for ii = 1:length(V)
                    %first, check that 1) there is something at index ii, 2) that
                    if(~isnan(P(ii).Centroid(1)) && A(ii) > minAreathresh && A(ii) < maxAreathresh && MI(ii) > minMIthresh && MI(ii) < maxMIthresh)
                        T_array(tt,1:2) = P(ii).Centroid; %store position in first two columns of T_array
                        T_array(tt,1:2) = T_array(tt,1:2) - [xshift, yshift]; %remove shift due to movement of the entire frame
                        T_array(tt,3)   = M(ii);          %store mass in third column
                        T_array(tt,4)   = time;           %store time from first frame in seconds in 4th column
                        T_array(tt,5)   = A(ii);          %store area in fifth column
                        T_array(tt,6)   = SF(ii);         %store shape factor in sixth column
                        T_array(tt,7)   = Rlistlab(ii);   %store RFP FL intensity in seventh column
                        T_array(tt,8)   = Glistlab(ii);   %store GFP FL intensity in eigth column
                        T_array(tt,9)   = listlabel(ii);  %store transferred mito intensity in ninth column
                        tt = tt+1;                        %increment T_array index
                    end
                end
            end
            clear Glistlab Rlistlab listlabel; %erase data from this step to prevent overwrite
        end
        
        %%
        
        if ~isempty(T_array) && sum(T_array(:,4) ~= T_array(1,4))>0
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Cell tracking using the track function
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            T_array(:,3) = T_array(:,3).*massfact; %change mass weighting in T_array
            T_array = sortrows(T_array, 4); %sort T_array based on time vectors
            minTx =  min(T_array(:,1));
            T_array(:,1) = T_array(:,1) -minTx +1; %make sure all x positions are positive\
            minTy =  min(T_array(:,2));
            T_array(:,2) = T_array(:,2) -minTy +1; %make sure all y positions are positive, new with rev6
            %move time to last column, T_array will now be [x, y, m, A, SF, t]
             T_array_temp = [T_array(:,1:3),T_array(:,7:9), T_array(:,5:6), T_array(:,4)];
            
            tracks = track(T_array_temp,max_disp,param);
            
            %move area back to 5th column, T_array will now be [x, y, m, t, A, SF] and
            %tracks will now be [x, y, m, t, cellnum, A, SF]
            tracks_temp = [tracks(:,1:3), tracks(:,9:10), tracks(:,7:8),tracks(:,4:6)];
            tracks = tracks_temp;
            
            T_array(:,3) = T_array(:,3)./massfact; %adjust mass weighting back to the way it was
            T_array(:,1) = T_array(:,1) + minTx -1; %set all x positions back to the way they were
            T_array(:,1) = T_array(:,2) + minTy -1; %set all y positions back to the way they were
            tracks(:,3) = tracks(:,3)./massfact; %adjust for mass weighting in the tracks array as well
            tracks(:,1) = tracks(:,1) +minTx -1; %set all x positions back to the way they were
            tracks(:,2) = tracks(:,2) +minTy -1; %set all y positions back to the way they were
        else
            tracks = [];
        end
        
        %save D_stored and L_stored in a separate .mat file to save memory
        
        %save(['data', row, col, '.mat'], 'D_stored', 'L_stored')
        parsave([savefolder 'data_RFP', num2str(Loc), '.mat'], D_stored, 'D_stored', 0)
        parsave([savefolder 'data_RFP', num2str(Loc), '.mat'], L_stored, 'L_stored', 1)
        
        %save tracks data and clear vector so that code can be parallelized
        parsave([savefolder 'Tdata_RFP', num2str(Loc), '.mat'], tracks, 'tracks', 0)
        parsave([savefolder 'Tdata_RFP', num2str(Loc), '.mat'], T_array, 'T_array', 1)
        parsave([savefolder 'Tdata_RFP', num2str(Loc), '.mat'], xshift_store, 'xshift_store', 1)
        parsave([savefolder 'Tdata_RFP', num2str(Loc), '.mat'], yshift_store, 'yshift_store', 1)
        parsave([savefolder 'Tdata_RFP', num2str(Loc), '.mat'], t_stored, 't_stored', 1)
        
    end %close if statement looking for stored data
end %close rows for loop

clear tracks T_array xshift_store yshift_store t_stored D_stored L_stored T_array_temp tracks_temp

% matlabpool close
%%
%reconstitute tracks vector
tracks_comp = [];
xshift_store_c = [];
yshift_store_c = [];
t_stored_c = [];
ii_stored = [];
maxindex = 0;
Loc_stored_c = [];
for Loc = 1:numLoc
    
    load([savefolder 'Tdata_RFP', num2str(Loc), '.mat']);
    if ~isempty(tracks)
        tracks(:,5) = tracks(:,5)+maxindex;
    end
    
    tracks_comp = [tracks_comp; tracks];
    xshift_store_c = [xshift_store_c; xshift_store(:,Loc)];
    yshift_store_c = [yshift_store_c; yshift_store(:,Loc)];
    t_stored_c = [t_stored_c; t_stored(:,Loc)];
    ii_stored = [ii_stored, 1:length(t_stored(:,Loc))'];
    Loc_stored_c = [Loc_stored_c, (1:length(t_stored(:,Loc))').*0 + Loc];
    
    if ~isempty(tracks_comp)
        maxindex = max(tracks_comp(:,5));
    end
    
    clear tracks xshift_store yshift_store t_stored
    
end

Loc_stored = Loc_stored_c;
tracks = tracks_comp;
xshift_store = xshift_store_c;
yshift_store= yshift_store_c;
t_stored = t_stored_c;
ii_stored = ii_stored';
clear tracks_comp xshift_store_c yshift_store_c t_stored_c maxindex Loc_stored_c

T0=min(tracks(:,4)); %find time of first image in the set
tracks(:,4) = (tracks(:,4)-T0);
t_stored = t_stored - T0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Code to extract growth rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FL signal specs paramater decided after completing tracking
% set these parameters based on intensity trends from tracking results
RFPInt=0.15;      % minimum RFP FL signal intensity to consider tumor cell with tagging
GFPInt=0.3;       % minimum GFP FL signal intensity to consider macrophage with emerald tagging
MitoInt=0.0165;   % minimum mito puncta intensity in tumor cell to count as puncta & not noise
MitoLng=4;        % minimum number of frames the puncta must exist in tumor cell to be counted as real
minpathlength=60; % logs tracks above minpathlength for growth calculation
%%
[num, indices] = track_numpart(tracks,minpathlength); %find tracks >= minpathlength
for ii = 1:num
    currentnum = tracks(indices(ii),5);
    [x,y,z,t] = track_partn_SF(tracks,currentnum); % calls x-y centroid coordinates(x&y), mass(z) and time(t)
    BF= polyfit(t, z, 1); %find best fit line to the data
    GRate(ii) = BF(1); %store growth rate
    Mass0(ii) = z(1); %store initial masss
    MassF(ii) = z(end); %store final mass
    Length(ii) = length(x); % stores length of tracks
end


%%
figure(2)
%plot growth over time
[num, indices] = track_numpart(tracks,minpathlength);
%plot those tracks

for ii = 1:num
    currentnum = tracks(indices(ii),5);
    [x,y,z,t,r,g,l] = track_partn_SF(tracks,currentnum);
    if mean(nonzeros(r))>=RFPInt && mean(nonzeros(g))<=GFPInt
       if sum(l>=MitoInt)>=MitoLng
           plot(t, z, '.-b')
           hold on
       else
           plot(t, z, '.-r')
           hold on
       end
    end
    if mean(nonzeros(g))>=GFPInt && mean(nonzeros(r))<=RFPInt
        plot(t, z, '.-g')
        hold on
    end
end
hold off

ylabel('Mass (pg)', 'FontSize', 14)
xlabel('time (h)', 'FontSize', 14)
set(gca, 'FontSize', 14)
title('labels: "cell #, initial mass, best fit slope, best fit intercept"', 'FontSize', 10)


%find and store data from this plot
h = findobj(gca,'Type','line');
fig2x=get(h,'Xdata');
fig2y=get(h,'Ydata');


%%
%save data
if ~exist([savefolder 'data_allframes_RFP']) || overwrite
    save([savefolder 'data_allframes_RFP'])
end

%% check one position tracking results to visually validate thge tracking results
load([savefolder 'data_RFP4'])  % random position selected
%%
figure(12)
for nn = 1:numf
D1F = imfilter(abs(D_stored(:,:,nn)), fspecial('gaussian', [10 10], 2));
% D1F = D_stored(:,:,nn);
% L2 = imagesegment_aggressive(D1F);
% [V1, M1, A1, MI1, P1, SF1] = imageprops_SF(L_stored(:,:,nn), D1F, pxlsize); %compute image properties
[V1, M1, A1, MI1, P1, SF1] = imageprops_SF(L_stored(:,:,nn), D1F, pxlsize); %compute image properties
%imagesc(D_stored(:,:,nn))
% plotBWoutlines(D_stored(:,:,nn),L_stored(:,:,nn),12);
plotBWoutlines(D_stored(:,:,nn),L_stored(:,:,nn),12);
colormap parula
title('areas in nth frame that are counted as "cells" labeled with masses')
lowA1 = find(A1 > minAreathresh & A1 < maxAreathresh & MI1 > minMIthresh & MI1 < maxMIthresh); %find indices of objects which are counted as "cells"
% lowA1 = find(A1); %uncomment this to select all the cells in the image
for k = 1 : length(lowA1)
    rectangle('Position', P1(lowA1(k)).BoundingBox, 'EdgeColor','y'); %plot yellow box around each cell
    text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(M1(lowA1(k))), 'Color', [0 1 0]) %label each cell with its mass in the first frame    text(P1(lowA1(k)).Centroid(1), P1(lowA1(k)).Centroid(2), num2str(M1(lowA1(k))), 'Color', [1 1 0]) %label each cell with its mass in the first frame
end
hold off;
pause(1); 
end


