%mean_cell_data_fitting_pred_fret_space.m written  9-25-17 by JTN to perform
%OLS optimization for the equation u_t + v(t)u_x =0 in a parallelized
%for loop to mean experimental data

%in this script, v(t) = v_1(m) + v_2(dm/dx) to look into chemokinesis and
%chemotaxis, where m(t,x) is interpolated FRET ratio data

function mean_cell_data_fitting_pred_fret_space(l,m,pred_ind,stat_model,simnum)

    %for bookkeeping
    well = m;
    
    %interp IC
    IC_type = 'interp';
    
        
    %select which model and which grid size we're using
    i = mod(l,4);
    if i == 0
        i = 4; %to determine grid size.
    end
    
    %load in data
    load('ind_cell_prof_data.mat')
    
    %choose data
    cell_data = avg_cell_data(ind_cell_data{m-1,2}(:,:,16:end),pred_ind)';
    FRET_data = avg_cell_data(ind_fret_data{m-1,2}(:,:,16:end),pred_ind)';
    m0 = 0;
    m1 = max(max(FRET_data));
    
        
    %initialize data grids
    [tndata,xndata] = size(cell_data);
    xdata = linspace(0,1,xndata);
    tdata = 5:1/3:1/3*(tndata-1+5*3);
    %make ndgrid for interpolation
    [Xdata,Tdata] = ndgrid(xdata,tdata);

    
    %generate grids for computation, also helpful for interpolation
    
    %grid sizes and models considered
    xnsize = [25 50 100 200];
    xn = xnsize(i);
    dt = 1e-3;
    [x,t] = grid_generate(xn,xdata(1),xdata(end),dt,tdata(1),tdata(end));
    x=x';
    tn = length(t);
    dx = x(2)-x(1);
    %interior, boundary points
    [x_int,xbd_0,xbd_1] = int_bd_def(xn);

    
    %%%%% Fit linear interpolant to (t,x) FRET data
    F1 = griddedInterpolant(Xdata,Tdata,FRET_data');
    % the vector x will be constant from now on, so can call interpolated
    % FRET data by: F(x,t*ones(length(x),1))
    
    %now, we can interpolate the derivative as well. First iterate through
    %each timepoint and smooth and take derivative:
    FRET_data_dx = zeros(size(FRET_data));
    for i = 1:tndata
        FRET_data_dx(i,:) = [diff(smooth(FRET_data(i,:)))' ...
            diff(FRET_data(i,end-1:end))]/(xdata(2)-xdata(1));
    end
    
    dm0 = min(min(min(FRET_data_dx)));
    dm1 = max(max(max(FRET_data_dx)));
    
    F2 = griddedInterpolant(Xdata,Tdata,FRET_data_dx');
    
    
    %%%% Now fit to migration data
    q_all = cell(simnum,1);
    q0_all = cell(simnum,1);
    J_all = zeros(simnum,1);

   

    %initial condition
        switch IC_type
        case 'interp'

                IC = interp1(xdata,cell_data(1,:),x);
                cutoff_x = leading_edge_calc(IC,x,.05,0);
                IC(x>cutoff_x)=0;

        case 'step'

                LE_loc = leading_edge_calc(smooth(cell_data(1,:)),xdata,0.8,0);
                IC = double(x<=LE_loc);
        end


    %boundary conditions
    BC_x_0 = @(t) 1;
    BC_x_1 = @(t) 0;

    %sparse matrix as a function for computation

    %A matrices now must be space-dependent
    A_pos = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],...
        [vw.*(-1+sw/2); (ve.*(1-1*se/2)-vw.*sw/2); (ve.*se/2)],xn,xn);
    
    A_pos_0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],...
        [ve.*(1-1*se/2); ve.*se/2],xn,xn);

    A_pos_1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[vw.*(-1+sw/2); ...
        (-vw.*sw/2)],total,total);



    A_neg = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],...
        [(-vw.*sw/2); (ve.*se/2+vw.*sw/2-vw); (ve-ve.*se/2)],xn,xn);

    A_neg_1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[(-vw.*sw/2); ...
        (vw.*sw/2-vw)],xn,xn);

    A_neg_0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],...
        [ve.*se/2; (ve-ve.*se/2)],xn,xn);

    options = optimset('maxiter',75);
    
    

    %each row of data matrix corresponds to data at given time point. In this
    %stat model, we only care about rows when t is a multiple of toCompare

    toCompare = 6;
    cell_data = cell_data(1:3*toCompare:end,:);

    tdata = tdata(1:3*toCompare:end);

    %get smaller data sample

    name_save = [stat_model '_' num2str(toCompare)];
        
    
    parfor k = 1:simnum

        
            if strcmp(IC_type,'interp') %don't estimate height
                %q = [v_1,v_2,..,v_8]^T;
                q0_all{k} = [.02*rand(4,1);.002*rand(4,1)];
                LB = [zeros(8,1)];
                UB = [inf(8,1)];
            else
                q0_all{k} = [.0075*rand(1,5),1];
                LB = [zeros(1,length(q0_all{k})-1),.7];
                UB = [inf*ones(1,length(q0_all{k})-1),1.2];
            end



            tic

            [q_all{k},J_all(k)] = fmincon(@(q) MLE_cost_D0_fewer_compare_space(cell_data...
                ,q,F1,F2,m0,m1,dm0,dm1,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,...
                tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg),q0_all{k}...
                ,[],[],[],[],LB,UB,[],options);

            toc

                       
    end

    save(['/scratch/summit/jona8898/FRET_fitting/FRET_interp_est_well_' ...
        num2str(well) '_' num2str(l) '_' name_save '_pred_'...
        num2str(pred_ind) '_space.mat' ],...
        'q_all','q0_all','J_all','p')

end
