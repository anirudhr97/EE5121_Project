function [y] = soft_thresholding(x,l)
    y = sign(x).*(max(abs(x)- l, 0));
end
