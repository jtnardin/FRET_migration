%FRET_dep_convection written 5-29-17 by JTN to simulate a 1d convection
%equation, where the rate of convection depends on the rate of MAPK
%activation

function udata = FRET_dep_convection_space(q_est,F1,F2,m0,m1,dm0,dm1,x,dx,xn,x_int,xbd_0,...
    xbd_1,t,dt,tn,tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg)

    %initialization, advection rates
   
     
    %FRET given by interpolation done previously
    FRET = @(t) F1(x,t*ones(xn,1));
    dFRETdx = @(t) F2(x,t*ones(xn,1));
    
    
    %points along m that we care about
    n=4;
    %for v(m)
    msamp1 = augknt([m0,m1,linspace(m0,m1,n)],2);
    %create spline functions
    v_spline1 = spmak(msamp1,q_est(1:4)');
    %for v(dm/dx)
    msamp2 = augknt([dm0,dm1,linspace(dm0,dm1,n)],2);
    %create spline functions
    v_spline2 = spmak(msamp2,q_est(5:end)');
    
    
    %evaluate V(m(t))
    Vx_c = @(t) (fnval(v_spline1,FRET(t)) + fnval(v_spline2,dFRETdx(t)))*dt/dx;

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

        
        v = Vx_c(t(i));
        
        %set BC
        u(xbd_0,i) = BC_x_0(t(i));
        u(xbd_1,i) = BC_x_1(t(i));

        %get Ax matrices
        [Axpos,Axneg] = AMatrixCompute(u(:,i-1),v,x_int,xn,A_pos,A_neg);

        u(:,i) = (speye(xn) + theta*(Axpos+Axneg))\(speye(xn) - ...
            (1-theta)*(Axpos+Axneg))*u(:,i-1);

        %[u(:,i),flag] = gmres((speye(xn*yn) + theta*(Ax + Ay)),(speye(xn*yn) - (1-theta)*(Ax + Ay))*u(:,i-1));


    end

    
    %interpolate model sim to match data
    
    [X,T] = meshgrid(x,t);
    
    udata = interp2(X,T,u',xdata,tdata');
        
end
