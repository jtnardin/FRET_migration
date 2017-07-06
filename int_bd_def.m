function [x_int,xbd_0,xbd_1] = int_bd_def(xn)

    xbd_0 = 1;
    xbd_1 = xn;

    xbd = union(xbd_0,xbd_1);
    
    x_int = 1:xn;
    
    x_int(xbd) = [];


end