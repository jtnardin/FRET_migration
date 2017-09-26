%positive_sensor.m written 9-16-16 by JTN to compute the sensors for flux
%limiters when the velocity is positive.

%updated 8-25-17 by JTN for the more general (2d) case. Currently not
%incorporating BC

% The code is a little messy because points corresponding to x (or y) =
% delta x (i.e., the secon row / column) will sample incorrectly, so we
% need to find these points and assign them a value of -1 (assuming 0
% Neumann BC)

% function [r_e,r_w,r_e_m0,r_w_m1] = positive_sensor(u,int,int_row_1,step)
function [r_e,r_w] = positive_sensor(u,int,int_row_1,step)

    %initialize r_w to have right number of points
    r_w = zeros(length(int),1);

    %compute r_e no trouble
    r_e = (u(int) - u(int-step))./(u(int+step) - u(int));    
    
    
    %only consider the points that will sample correctly (which is anything
    %not in the first row)
    int_past_row_1 = ~int_row_1;

    
    %For r_w, compute most points as normal
    r_w(int_past_row_1) = (u(int(int_past_row_1)-step) - u(int(int_past_row_1)-2*step))./...
        (u(int(int_past_row_1)) - u(int(int_past_row_1)-step));
        
    %and now fix points that would have sampled incorrectly
    %we are currently assuming a zero neumann BC
    r_w(int_row_1) = -1;
    
%   
% %     boundary
%     r_w_m1 = r_e(bd_nm1);
%     r_e_m0 = r_w(bd_0);
    
    
    %eliminate NaN values (0/0 -- not steep!)
    r_e(isnan(r_e)) = 1;
    r_w(isnan(r_w)) = 1;
%     r_w_m1(isnan(r_w_m1)) = 1;
%     r_e_m0(isnan(r_e_m0)) = 1;
%     
    %set inf values to large value
    r_e(isinf(r_e)) = 100;
    r_w(isinf(r_w)) = 100;
%     r_w_m1(isinf(r_w_m1)) = 100;
%     r_e_m0(isinf(r_e_m0)) = 100;
    
end