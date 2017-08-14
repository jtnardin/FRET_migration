clear all; clc

for well = 4:5



    sim = 3;
    pred_ind = 1;

    load(['FRET_interp_est_well_' num2str(well) '_fewer_6_pred_'...
        num2str(pred_ind) '_final.mat'])

    save_im = 1;
    plot_res = 1;
    stat_model = 'fewer';


    %options for each type of simulation

    run_sim(sim,well,pred_ind,q_final{sim},save_im,plot_res,stat_model)

%     pred_data(sim,well,pred_ind,q_final{sim},save_im,plot_res,stat_model)

end