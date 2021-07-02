function [h] = plotBWoutlines(I, L, fignum)
%[h] = plotBWoutlines(I, L)
%plots 
%returns figure handle

if nargin > 2
    if fignum~=-1
        h = figure(fignum);
    end
else
   h = figure;
end

BWoutline = bwperim(L>=1);
BWedge = bwboundaries(L>=1);
imagesc(I)
colormap gray
hold on
for k=1:length(BWedge)
    boundary = BWedge{k};
    plot(boundary(:,2), boundary(:,1),'r','LineWidth',2);
end
hold off