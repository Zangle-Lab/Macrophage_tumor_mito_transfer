function [ D, L,taglist,listGFP ] = LoadSegment_GFP_AlFL( fname, wavelength,GFP, limG)
%function [ D, L, B ] = LoadSegment( fname, wavelength )
%function to load and segment the data stored in fname
%phase data should be stored in fname.Phase

%step 1, load phase images, segment and label cells
Loaded = load(fname);
D = -Loaded.Phase;
D=D(139:486,62:512);                    %crop phase image to matc FL image area
B = imagebackground_poly4(D);           %compute image background
D = ((D-B).*wavelength);                %subtract background from data files
DF = imfilter(abs(D), fspecial('gaussian', [3 3], 2)); % blur phase image to prevent over-segmentation
[L,~] = imagesegment_aggressive(DF);    %segment image (detect distinct cell regions and disconnect connected regions)

%step 2, prepare separate masks for macrophage and tumor transferred mito FL signals
GFP = imfilter(GFP(:,1:669), fspecial('gaussian', [5 5], 1)); %crop FL image to overlap phase image area
GFPmask1=GFP>limG(1);
GFPmask2=GFP<limG(2);
GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*limG(2))+((1-GFPmask1).*limG(1));
GFP2=(GFP1-limG(1))/(limG(2)-limG(1));   % normalize GFP image
% rolling ball filter to separate high sptial frequency mito punctas from
% intense spread out macrophage FL signal:
GFP3=imtophat(GFP2, strel('sphere',8));  % change size of sphere based on average size of mito punctas
GFP4=imresize(GFP2,[348,451]);           % mask for macrophage FL signal
GFP5=imresize(GFP3,[348,451]);           % mask for transferred mito puncta FL signal

%step 3, log macrophage & transferred mito FL intensities
listLo=unique(single(L));
if sum(listLo)~=0
    for uu=2:length(listLo)
        taglist(uu-1,1)=listLo(uu);      % list of macrophage FL intensity
        taglist(uu-1,2)=mean(mean(nonzeros((L==listLo(uu,1)).*GFP4)));
    end
    for labn=2:length(listLo)
        listGFP(labn-1,2)=mean2(nonzeros((L==listLo(labn,1)).*GFP5));
        listGFP(labn-1,1)=listLo(labn);  % list of mito puncta FL intensity
    end
else
    taglist=[0,0]; listGFP=[0,0];        % in case there is no cell in frame
end

end