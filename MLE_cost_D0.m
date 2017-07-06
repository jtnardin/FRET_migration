function [J,WLS_SV,weight_f,weight_matrix,res,model] = MLE_cost_D0(cell_data,q_est,p,m0,m1,...
    x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg,weightno)

    weight_f = zeros(weightno,1);
    N = numel(cell_data);
    
    [model] = FRET_dep_convection(q_est,p,m0,m1,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,IC,IC_type,...
        BC_x_0,BC_x_1,A_pos,A_neg);
    
    %keeps track of which data entries correspond to which weights
    %(1,2,...,11)
    model_subset_ind = zeros(numel(cell_data),1);
    %a vector full of the weights that correspond to points
    weight_matrix = zeros(numel(cell_data),1);
    
    %calculate residuals
    res = model (:) - cell_data(:);
    
    OLS_SV = 1/N*sum(res.^2); %slightly biased sample variance
    
    %iterate through number of weights needed
    for j = 1:weightno

                %different weight numbers need to select from different
                %areas of the data.
                if weightno == 11
                                %where is model between different densities
                    model_subset = (model(:) >= .1*(j-1))&(model(:) < .1*(j)); 

                elseif weightno == 29

                    if j <= 10 %10 bins between 0 and 0.1
                        model_subset = (model(:) >= .01*(j-1))&(model(:) < .01*(j));
                    elseif j >=11 && j <= 18 %8 bins between 0.1 and 0.9
                        model_subset = (model(:) >= .1*(j-10))&(model(:) < .1*(j-9));
                    elseif j >= 19 && j <= 28 %10 bins between 0.9 and 1
                        model_subset = (model(:) >= .9 + (j-19)*.01)&(model(:) < .9 + (j-18)*.01);
                    elseif j == 29 %in case data ever above 0.1
                        model_subset = (model(:) >= 1);
                    end

                end

                %take corresponding subset of res
                res_subset = res(model_subset); 
                %calculate weight for that section
                weight_f(j) = sqrt((1/numel(res_subset)*sum(res_subset.^2))/OLS_SV); 

                %keep track of which points correspond to which weights
                model_subset_ind(model_subset) = j;
                %put weights into corresponding location in matrix
                weight_matrix(model_subset) = weight_f(j);
    end


 
   
    WLS_SV = 1/N*sum((res./weight_matrix).^2);
    
    J = sum(log(sqrt(WLS_SV)*weight_matrix) + (res./weight_matrix).^2/(2*WLS_SV));


end