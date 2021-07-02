function [M] = backgroundmask(I)
%function [M] = backgroundmask(I)
%function to mask the background of an image, I
%input: I, the grayscale image to find the background of
%output: B, the background of I
%method: find 'objects' in I, mask them from the image, fit the remaining
%pixels with an ellipse, subtract that from the image

[junk threshold] = edge(I, 'sobel');
fudgeFactor = .4;
%fudgeFactor = 0.8;
BWs = edge(I,'sobel', threshold * fudgeFactor);

%step 2, dilate the image
se90 = strel('line', 8, 90);
se0 = strel('line', 8, 0);
BWsdil = imdilate(BWs, [se90 se0]);
BWdfill = imfill(BWsdil, 'holes');

BWf = imerode(BWdfill, strel('disk', 2));

M = bwareaopen(BWf,50); %remove all areas of fewer than 50 pixels


end



