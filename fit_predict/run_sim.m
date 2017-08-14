function run_sim(l,m,pred_ind,q_est,save_im,res_plot,stat_model)


    i = mod(l,4);
    if i == 0
        i = 4; %to determine grid size.
    end
    
    %load in data
    load('ind_cell_prof_data.mat')
    load('FRET_data_ind.mat')

    %choose data
    cell_data = avg_cell_data(ind_cell_data{m-1,2}(:,:,16:end),pred_ind)';
    FRET_data = avg_fret_data(FRET_time_ind{m-1}(:,16:end),pred_ind);
    
    %initialize data grids
    [tndata,xndata] = size(cell_data);
    xdata = linspace(0,1,xndata);
    tdata = 0:1/3:1/3*(tndata-1);

    %IC type
    IC_type = 'interp';
    weightno = 29;
    
    
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


    tic
    
%     [model] = FRET_dep_convection(q_est,p,m0,m1,x,dx,xn,x_int,xbd_0,xbd_1,t,...
%         dt,tn,tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg);
%     

    if strcmp(stat_model,'WLS')

        [J,WLS_SV,weight_f,weight_matrix,res,model] = MLE_cost_D0_WLS(cell_data,q_est,p,m0,m1,...
            x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg,weightno);
        
    elseif strcmp(stat_model,'fewer')
        
        %each row of data matrix corresponds to data at given time point. In this
        %stat model, we only care about rows when t is a multiple of toCompare
                
        toCompare = 6;
        cell_data = cell_data(1:3*toCompare:end,:);
        
        tdata = tdata(1:3*toCompare:end);
        
        %get smaller data sample
        
        name_save = [stat_model '_' num2str(toCompare)];
        
        [J,OLS_SV,res,model] = MLE_cost_D0_fewer_compare(cell_data,q_est,p,m0,m1,...
            x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg);
        
        
    elseif strcmp(stat_model,'autor')
        
        
        %size of data
        xnd = length(xdata);
        tnd = length(tdata);

        
        %set alpha
        alpha = .4;
        alpha_vec = [repmat([alpha*ones(1,xnd-1) 0],1,tnd-2) alpha*ones(1,xnd-1)];
        
        total = xnd*tnd;
        off_diag = xnd+2:total;
          
        A = sparse([1:total off_diag],[1:total off_diag-xnd-1],[ones(1,total) alpha_vec],total,total); 
        stat_options.AinvAT = inv(A*A');
        stat_options.Ainv = inv(A);
        
        [J,res,model] = MLE_cost_D0_autor(cell_data,q_est,p,m0,m1,...
            x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg,stat_options);
        
        name_save = [stat_model '_alpha' num2str(alpha)];
    
    end


    toc
    
    figure
    hold on
    
    colors = distinguishable_colors(length(tdata));
    
    if strcmp(stat_model,'fewer')
        
        for i = 1:length(tdata)
            plot(xdata,model(i,:),'color',colors(i,:))
            plot(xdata,cell_data(i,:),'.','color',[colors(i,:)])
        end
        
    else
    
        count = 1;
        
        for i = 1:floor(tndata/4):tndata
            plot(xdata,model(i,:),colors(count))
            plot(xdata,cell_data(i,:),[colors(count) '.'])

            count = count + 1;
        end
    end
    
%     legend('1','2','3','4','5','location','northeast')
    
    cell_dens = [1700,2500,3000,4000];

    h=title(['Fitting cell profile, density = ' num2str(cell_dens(m-1)) ' cells/mm$^2, x_n$ = ' num2str(xn)]);
    set(h,'interpreter','latex')
    xlabel('Location (x)')
    ylabel('u(t,x)')
    
    if save_im == 1
        exportfig(gcf,['cell_fitting_' num2str(m) '_' num2str(l) name_save '_pred_' num2str(pred_ind) '.eps'],'color','rgb')
        saveas(gcf,['cell_fitting_' num2str(m) '_' num2str(l) name_save '_pred_' num2str(pred_ind) '.fig'])
    end
    
%     close
    
%     figure
%     hold on
%     
%     plot(tdata,ppval(p,tdata),'b')
%     
%     if plot_data == 1
%         plot(tdata,FRET_data(1:tndata),'b.')
%         
%     end
%     
%     title(['FRET interpolation for cell density = ' num2str(cell_dens(m-1)) ' cells/mm$^2$'],'interpreter','latex')
%     xlabel('Time (t)')
%     ylabel('FRET')
%     
%     if save_im == 1
%         exportfig(gcf,['FRET_interp_' num2str(m) '_' name save '.eps'],'color','rgb')
%         saveas(gcf,['FRET_interp_' num2str(m) '_' name save '.fig'])
%     end
    
%     close

    if res_plot == 1
       
       [X,T] = meshgrid(xdata,tdata);
       
       
       if strcmp(stat_model,'WLS')
            res_final = (res./weight_matrix)/sqrt(2*WLS_SV);
       elseif strcmp(stat_model,'autor')
           res_final = stat_options.Ainv*res/sqrt(N);
%             res_final = reshape(res_final,)
           figure

           surf(X,T,reshape(res_final,length(tdata),length(xdata)),'edgecolor','none')

           xlabel('x','fontsize',30)
           ylabel('t','fontsize',30)

           set(gca,'fontsize',20)

           title(['WLS residuals, cell density = ' num2str(cell_dens(m-1)) ' cells/mm$^2$'],'interpreter','latex')
           axis([0 1 0 tdata(end)])

           caxis([-2 2])
           colorbar
       
       
            if save_im == 1
                exportfig(gcf,['FRET_interp_' num2str(m)   '_' num2str(l) '_' name save '_res.eps'],'color','rgb')
                saveas(gcf,['FRET_interp_' num2str(m) '_' num2str(l) '_' name save '_res.fig'])
            end
    
       elseif strcmp(stat_model,'fewer')
           
           model_vec = model(:);
           
           res_final = res;

           figure
           hold on
           plot(res_final,'b.')
           plot([0 1e6],[0 0],'k')
           axis([0 length(model_vec) -.3 .3])

           xlabel('Data Points','fontsize',30)
           ylabel('$y_{ij} - f(t_i,x_j)$','fontsize',30,'interpreter','latex')
           title(['Fit residual, cell density = ' num2str(cell_dens(m-1)) ' cells/mm$^2,\ x_n$ = ' num2str(xn)],'interpreter','latex')


           if save_im == 1
             exportfig(gcf,['FRET_interp_fit_' num2str(m)   '_' num2str(l) '_' name_save '_res_pred_' num2str(pred_ind) '.eps'],'color','rgb')
             saveas(gcf,['FRET_interp_fit_' num2str(m) '_' num2str(l) '_' name_save '_res_pred_' num2str(pred_ind) '.fig'])
           end
           
       end
        
    end
    
end