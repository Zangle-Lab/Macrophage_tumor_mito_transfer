function [ishift, jshift] = CorrShift(A,B)
%find correlation shift between two matrices of the same size, A and B
%[ishift,jshift] is the amount to move elements of A such that they line up
%with elements of B
C = normxcorr2(A,B);

idim = length(A(:,1));
jdim = length(A(1,:));

[YY,II] = max(C);
[Y,J] = max(YY);

ishift = II(J)-idim;
jshift = J-jdim;