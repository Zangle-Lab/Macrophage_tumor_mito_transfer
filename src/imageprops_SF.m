function [V, M, A, MI, P, SF] = imageprops_SF(L, I, pxlsize)
%function [V, M, A, MI, P, SF] = imageprops_SF(L, I, pxlsize)
%function to return volume, mass and area of regions L in image I
%inputs: L, the label image *can also be BW mask image. should be nonzero
%in regions where the image will be analyzed; I, the image to be processed;
%pxlsize, the size of each pixel in the image (mm/pixel)
%outputs: V, the measured volume (um^3) of each region; M the measured mass
%(pg); A, the area in pixels of each region; MI, the mean intensity of each
%region; P, struct containing all data returned by the call to the
%regionprops function ('Area', 'Centroid', 'MeanIntensity', 'BoundingBox',
%'PixelValues') and the standard deviation ('Std')
%_SF version also returns 'shape factor' defined as the circularity
%(= 4pi*A/P^2)


%%% define relationship between mass and "volume"
K = 1./(10000).^3./100./0.0018.*1e12; %pg/um^3

P = regionprops(L, I, 'Area', 'Centroid', 'MeanIntensity', 'BoundingBox', 'PixelValues');
Mask = backgroundmask(I); %find mask with more restrictive areas
L2 = L;
L2(~Mask) = 0;
L2 = imkeeplargest(L2);
P2 = regionprops(L2, 'Area', 'Perimeter');
% yy=0;
if max(max(L))>0
    for ii = 1:length(P)
%         if P(ii).Area>0
%             yy=yy+1;
            yy=ii;
            A(yy) = P(ii).Area;
            MI(yy) = P(ii).MeanIntensity;
%             Ct(yy,1:2)= P(ii).Centroid;
            P(ii).Std = std(P(ii).PixelValues);
            if ii<=length(P2) && ~isempty(P2(ii))
                SF(yy) = P2(ii).Area*4*pi/(P2(ii).Perimeter.^2);
            else
                SF(yy) = NaN;
            end
%         end
    end
else
    A = [];
    MI = [];
    P = [];
    SF = [];
end

V = MI.*A.*pxlsize.^2.*1e3;
M = V.*K;