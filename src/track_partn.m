function [x,y,z,t,r,g,l] = track_partn(tracks, n)
%function [x,y,z,t] = track_partn(tracks, n)
%returns the coordinates and times for particle n in the tracks array
%generated by the track.m function
%inputs: tracks, array containing x, y, z, t, i (locations x, y, third
%coordinate z, and time coordinate, t for i particles); n, particle number
%to return
%outputs: x, y, z, t: coordinates of particle n
%TAZ 4/27/10

if n > max(tracks(:,5))
    disp('invalid particle number')
    x = [];
    y = [];
    z = [];
    t = [];
    r = [];
    g = [];
    l = [];
else
    x = tracks(tracks(:,5) == n, 1);
    y = tracks(tracks(:,5) == n, 2);
    z = tracks(tracks(:,5) == n, 3);
    t = tracks(tracks(:,5) == n, 4);
    r = tracks(tracks(:,5) == n, 8);
    g = tracks(tracks(:,5) == n, 9);
    l = tracks(tracks(:,5) == n, 10);
end