%MLE_cost_D0_autor written 7-13-17 by JTN to do 29-point WLS scheme

%%7-13-17 still need formula for sigma^2

function [J,res,model] = MLE_cost_D0_autor(cell_data,q_est,p,m0,m1,...
    x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg,stat_options)

    AinvAT = stat_options.AinvAT;
   
    %run simulation
    [model] = FRET_dep_convection(q_est,p,m0,m1,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,IC,IC_type,...
        BC_x_0,BC_x_1,A_pos,A_neg);
    
    
    %each row of model corresponds to solution at given time point. In this
    %stat model, we only care about rows when t is a multiple of toCompare
    %same with data;
    
%     model = model(mod(tdata,toCompare)==0,:);
%     cell_data = cell_data(mod(tdata,toCompare)==0,:);
    
    %total data points considered
    N = numel(cell_data);
    
    model = reshape(model',N,1);
    cell_data = reshape(cell_data',N,1);
    
    %calculate residuals
    res = model(:) - cell_data(:);
    
    
    %%%%%need to update and calculate SV
%     OLS_SV = 1/N*sum(res.^2); %slightly biased sample variance
    
    
    J=1/N*res'*AinvAT*res;


end