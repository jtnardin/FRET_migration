%create x, t grids for analysis

function [xgrid,tgrid] = grid_generate(xn,x0,xend,dt,t0,tend)


    xgrid = linspace(x0,xend,xn);

    tgrid = t0:dt:tend+dt;


end