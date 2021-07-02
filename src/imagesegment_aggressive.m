function [L,BWdfill] = imagesegment_aggressive(I)
%function to segment an image of cells and label them
%procedure followed by breaking up connected cells using a watershed image transform
%input: I, the grayscale image to segment
%output: L, the labeled, segmented image

% step 1, detect cells by edge features
[junk threshold] = edge(I, 'sobel');
fudgeFactor = 0.4; %was 0.4 for RBCs
BWs = edge(I,'sobel', threshold * fudgeFactor);

%step 2, dilate the image
se90 = strel('line', 10, 90);
se0 = strel('line', 10, 0);
BWsdil = imdilate(BWs, [se90 se0]);

%step 3, fill gaps
BWdfill = imfill(BWsdil, 'holes');

%step 4, smooth image
seD = strel('diamond',2);
BWfinal = imerode(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 500); %remove regions smaller than 500 pixels, debris

% step 5, disconnecting cells by watershed
Igr = mat2gray(I);
thr = 0.11; % threshold to determine valleys (cell high areas points) for watershed, varies with phase data & cell type
mask_em = imextendedmax(Igr, thr);
%Next step: complement the image so that the peaks become valleys. We do
%this because we are about to apply the watershed transform, which
%identifies low points, not high points. 
I_c = imcomplement(I);
%Next: modify the image so that the background pixels and the extended
%maxima pixels are forced to be the only local minima in the image. 
I_mod = imimposemin(I_c, ~BWfinal | mask_em);
%now, compute watershed transform
L = watershed(I_mod);

%step 6, remove connected images on border
L = imclearborder(L, 4);
