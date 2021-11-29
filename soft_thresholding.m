function [y] = soft_thresholding(x,l)
%{
Function to implement the soft thresholding used in the ISTA and FISTA.
%}
    y = sign(x).*(max(abs(x)- l, 0));
end
