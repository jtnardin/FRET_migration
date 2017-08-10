clear all; clc

well = 2;

sim = 3;
pred_ind = 1;

load(['FRET_interp_est_well_' num2str(well) '_fewer_12_pred_'...
    num2str(pred_ind) '_final.mat'])

save_im = 0;
plot_res = 1;
stat_model = 'fewer';


%options for each type of simulation


pred_data(sim,well,pred_ind,q_final{sim},save_im,plot_res,stat_model)

    