%diff_conv.m written 1-26-17 by JTN to run a 1d diffusion-advection
%equation

function [udata,mdata] = diff_conv_f_vel_spline_OLS_D0(q_est,q_calc,x,dx,xn,x_int,xbd_0,...
    xbd_1,t,dt,tn,tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg)

    %initialization, advection rates
   
%     q_est = [alpha1,alpha2,beta,gamma]
%     q_calc = [a,b];
     
    FRET = @(t) ept_fcn_pulse(q_est,q_calc,t);
        
    Vx_c = @(t) q_est(1)*FRET(t)*dt/dx;


    %crank nicholson
    theta = 0.5;

   
    %sigma for flux limiters
    sigma = @(r) (r+abs(r))./(1+abs(r));


    %initialize
    u = zeros(xn,tn);
   
    %initial condition --     
    if strcmp('step',IC_type) %step function with variable height
    %drops at calculated leading edge.
        u(:,1) = q(end)*IC;
    elseif strcmp('interp',IC_type) %interp first
        %data point
        u(:,1) = IC;
    end


    for i = 2:tn


        %set BC
        u(xbd_0,i) = BC_x_0(t(i));
        u(xbd_1,i) = BC_x_1(t(i));

        %get Ax matrix

        if Vx_c(t(i)) >= 0

            r_e = (u(x_int,i-1) - u(x_int-1,i-1))./(u(x_int+1,i-1) - u(x_int,i-1));    
            %for r_w, start at yn-1 to avoid sampling off the grid ... assign
            %value of -1 for first values (would be true for 0 Neumann BC)
            r_w = (u(x_int(2:end)-1,i-1) - u(x_int(2:end)-2,i-1))./(u(x_int(2:end),i-1) - u(x_int(2:end)-1,i-1));
            r_w = [-1;r_w];



            %eliminate NaN values (0/0 -- not steep!)
            r_e(isnan(r_e)) = 1;
            r_w(isnan(r_w)) = 1;
            %r_w_m1(isnan(r_w_m1)) = 1;
            %r_e_m0(isnan(r_e_m0)) = 1;

            %set inf values to large value
            r_e(isinf(r_e)) = 100;
            r_w(isinf(r_w)) = 100;


            Ax_exp = A_pos(sigma(r_e),sigma(r_w),Vx_c(t(i-1)),x_int,1);
            Ax_imp = A_pos(sigma(r_e),sigma(r_w),Vx_c(t(i)),x_int,1);


        elseif Vx_c(t(i)) < 0

            %remove last x strip to avoid sampling off grid ... replace with -1
            r_e = (u(x_int(1:end-2)+1,i-1) - u(x_int(1:end-2)+2,i-1))./(u(x_int(1:end-2),i-1) - u(x_int(1:end-2)+1,i-1));
            r_e = [r_e;-1];
            r_w = (u(x_int,i-1) - u(x_int+1,i-1))./(u(x_int-1,i-1) - u(x_int,i-1));


            %eliminate NaN values (0/0 -- not steep!)
            r_e(isnan(r_e)) = 1;
            r_w(isnan(r_w)) = 1;
            %r_w_m1(isnan(r_w_m1)) = 1;
            %r_e_m0(isnan(r_e_m0)) = 1;

            %set inf values to large value
            r_e(isinf(r_e)) = 100;
            r_w(isinf(r_w)) = 100;

            Ax = A_neg(sigma(r_e),sigma(r_w),Vx_c,xy_int,yn);

        end



        u(:,i) = (speye(xn) + theta*Ax_imp)\(speye(xn) - (1-theta)*Ax_exp)*u(:,i-1);

        %[u(:,i),flag] = gmres((speye(xn*yn) + theta*(Ax + Ay)),(speye(xn*yn) - (1-theta)*(Ax + Ay))*u(:,i-1));


    end

    
    %interpolate model sim to match data
    
    [X,T] = meshgrid(x,t);
    
    udata = interp2(X,T,u',xdata,tdata');
    mdata = ept_fcn_pulse(tdata);
    
end
