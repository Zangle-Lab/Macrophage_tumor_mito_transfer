function [ D, L,taglist] = LoadSegment_RFP_AlFL( fname, wavelength,RFP,limR)

% section 1: segment to make phase image cell label
Loaded = load(fname);                   % phase image filename is fname 
D = -Loaded.Phase;
D=D(139:486,62:512);                    % cropping phase image to match FL image area
B = imagebackground_poly4(D);           %compute image background
D = ((D-B).*wavelength);                %subtract background from data files
DF = imfilter(abs(D), fspecial('gaussian', [3 3], 2)); % slightly blurring phase image reduces over-segmentation of cells 
L = imagesegment_aggressive(DF);        %segment image (detect distinct cell regions and disconnect connected regions)

% section 2: creating mask for FL image to log tumor cell RFP intensities
% FL image is cropped to overlap cropped phase image area
RFP = imfilter(RFP(:,1:669), fspecial('gaussian', [5 5], 1)); 
FPmask1=RFP>limR(1);
FPmask2=RFP<limR(2);
FP1=(single(RFP).*FPmask1.*FPmask2)+((1-FPmask2).*limR(2))+((1-FPmask1).*limR(1));
FP2=(FP1-limR(1))/(limR(2)-limR(1));    % normalize FL image
% FL image is higher pixel resolution, so resized to match phase image resolution
BWsn=imresize(FP2,[348,451]);
listL=unique(single(L));                % listing cell labels
if sum(listL)~=0
    for uu=2:length(listL)
        taglist(uu-1,1)=listL(uu);      % logging red Fl intensity of each tumor cell
        taglist(uu-1,2)=mean(mean(nonzeros((L==listL(uu,1)).*BWsn)));
    end
else
    taglist=[0,0];                      % in case no cells in frame
end

end