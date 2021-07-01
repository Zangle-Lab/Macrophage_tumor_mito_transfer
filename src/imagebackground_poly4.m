function [B,M] = imagebackground_poly4(I)
%function to find the background of an image, I, using a 6th order
%polynomial fit to the background pixels
%input: I, the grayscale image to find the background of
%output: B, the background of I
%method: find 'objects' in I, mask them from the image, fit the remaining
%area using polynomial fitting

%step 2, detect cells by edge feature
[junk threshold] = edge(I, 'sobel');
fudgeFactor = 0.7; %was 0.4 for RBCs
BWs = edge(I,'sobel', threshold * fudgeFactor);

%step 3, dilate the image
se90 = strel('line', 6, 90);
se0 = strel('line', 6, 0);
BWsdil = imdilate(BWs, [se90 se0]);

%step 4, fill gaps
BWdfill = imfill(BWsdil, 'holes');

%step 5, smooth image
seD = strel('diamond',1);
BWfinal = imerode(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 200); %remove regions smaller than 200 pixels (debris, not cell)

%step 6, polynomial fitting the regions except cells
IList = I(~BWfinal);

sz = size(I);
[X,Y] = meshgrid(1:sz(2), 1:sz(1));
XList = X(~BWfinal);
YList = Y(~BWfinal);

if sum(~isnan(BWfinal)) ~=0
    CFit = polyfitn([XList,YList], IList, 6);
    B =(reshape(polyvaln(CFit, [X(:), Y(:)]), sz(1), sz(2)));
else
    B = I;
end

M = ~BWfinal; %return the mask used for processing