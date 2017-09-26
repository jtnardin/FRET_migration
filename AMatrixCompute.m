%AMatrixCompute written 9-26-17 by JTN to create the negative and positive
%matrices for computation, given the velocity.

function [Axpos,Axneg] = AMatrixCompute(u,v,x_int,xn,A_pos,A_neg)

    
    
    ve_int = v(x_int+1);
    vw_int = v(x_int-1);
    

    %estimate positive, negative velocity locations
    V_pos_loc = ve_int + vw_int >= 0;
    V_neg_loc = ve_int + vw_int <  0;


    %compute velocities at neg, pos points
    vw_pos = vw_int(V_pos_loc);
    ve_pos = ve_int(V_pos_loc);

    vw_neg = vw_int(V_neg_loc);
    ve_neg = ve_int(V_neg_loc);
    
    

    %sigma for flux limiters
    sigma = @(r) (r+abs(r))./(1+abs(r));


    
    %locate x_int points for positive x vel
    x_intp = x_int(V_pos_loc); 
    %which points are on the second row?
    x_intp_row1 = mod(x_intp,xn)==2;
    %locate x_int points for negative x vel
    x_intn = x_int(V_neg_loc);
    %which points are in penultimate row?
    x_intn_rown = mod(x_intn,xn)==(xn-1);
   
    %create sensors
    [u_ep,u_wp] = positive_sensor(u,x_intp,x_intp_row1,1);
    [u_en,u_wn] = negative_sensor(u,x_intn,x_intn_rown,1);
        
    
    
    %compute Axpos, Axneg
    Axpos = A_pos(sigma(u_ep),sigma(u_wp),ve_pos,vw_pos,x_intp,1);
    Axneg = A_neg(sigma(u_en),sigma(u_wn),ve_neg,vw_neg,x_intn,1);
    
    
end