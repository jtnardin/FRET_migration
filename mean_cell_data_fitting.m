%mean_cell_data_fitting.m written  4-5-17 by JTN to perform WLS
%MLE optimization for the equation u_t = Du_xx - v(t)u_x in a parallelized
%for loop to mean experimental data

function mean_cell_data_fitting(l,m)

    %for bookkeeping
    welllet = 'mean';
    well = m;
    
    %interp IC
    IC_type = 'interp';
    %how many weights for WLS
    weightno = 29;
    
    %how many rounds of optimization do we want to do?
    simnum = 10;

    %select with model and which grid size we're using
    
    
    i = mod(l,4);
    if i == 0
        i = 4; %to determine grid size.
    end
    
    %load in data
    load('mean_cell_prof_data.mat')
    load('FRET_WT_mean_over_time.mat')

    %choose data
    cell_data = mean_cell_data{m-1,2}';
    FRET_data = FRET_mean_dens{m-1};
    FRET_data = FRET_data/max(FRET_data); %normalize
    
    %initialize data grids
    [tndata,xndata] = size(cell_data);
    xdata = linspace(0,1,xndata);
    tdata = 0:1/3:1/3*(tndata-1);

    
    %%%%% Fit hermite interpolant to FRET data
    %how much of the data to take (gets too jagged if I include all data
    %points, N = 10 works well for wells 2 and 3. Slightly off for wells 4
    %and 5, so going with N=5).
    N=5;

    FRET_data = FRET_data(1:tndata);
    p = pchip(tdata(1:N:end),FRET_data(1:N:end));
    m0 = min(FRET_data);
    m1 = max(FRET_data);
    %can evaluate p over time as ppval(p,tq), where tq are query points
    
    
    
    %%%% Now fit to migration data
       
    %grid sizes and models considered
    xnsize = [25 50 100 200];
    
    q_all = cell(simnum,1);
    q0_all = cell(simnum,1);
    J_all = zeros(simnum,1);

   

    %generate grids for computation
    xn = xnsize(i);
    dt = 1e-3;

    [x,t] = grid_generate(xn,xdata(1),xdata(end),dt,tdata(1),tdata(end));
    tn = length(t);
    dx = x(2)-x(1);
    %interior, boundary points
    [x_int,xbd_0,xbd_1] = int_bd_def(xn);

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

    A_pos = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-v+v.*sw/2); ...
        (v-v.*se/2-v.*sw/2); (v.*se/2)],xn,xn);

    A_neg = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-v.*sw/2); ...
        (v.*se/2+v.*sw/2-v); (v-v.*se/2)],xn,xn);

    options = optimset('maxiter',75);

    parfor k = 1:simnum

        
            if strcmp(IC_type,'interp') %don't estimate height
                %q = [v_1,v_2,..,v_5]^T;
                q0_all{k} = [.002].*rand(5,1);
                LB = zeros(5,1);
                UB = inf(5,1);
            else
                q0_all{k} = [.0075*rand(1,5),1];
                LB = [zeros(1,length(q0_all{k})-1),.7];
                UB = [inf*ones(1,length(q0_all{k})-1),1.2];
            end



            tic

            [q_all{k},J_all(k)] = fmincon(@(q) MLE_cost_D0(cell_data,q,p,m0,m1,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,...
                         IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg,weightno),q0_all{k},[],[],[],[],LB,UB,[],options);

            toc

                       
    end

    save(['/scratch/summit/jona8898/FRET_fitting/FRET_interp_est_' welllet num2str(well) '_' num2str(l) '.mat' ],...
                'q_all','q_fret','q0_all','J_all')

end
