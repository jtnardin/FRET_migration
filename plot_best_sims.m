 clc

for well = 2:5
    for sim = 1:3


        load(['FRET_interp_est_mean' num2str(well) '_final.mat'])


        plot_data = 1;
        save_im = 1;
        plot_res = 1;
        
        run_sim(sim,well,q_final{sim},plot_data,save_im,plot_res)
    end
end