%FRET_interpolation_1d.m written 7-6-17 by JTN to do a spline or hermite
%interpolation of average FRET data over time

%which well?
m = 4;

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

plot(tdata,FRET_data(1:length(tdata)))

t_fine = linspace(0,tdata(end),200);

N=5;

p = pchip(tdata(1:N:end),FRET_data(1:N:length(tdata)));
s = spline(tdata(1:N:end),FRET_data(1:N:length(tdata)));

plot(tdata,FRET_data(1:length(tdata)))
hold on
plot(t_fine,ppval(p,t_fine),'k')
plot(t_fine,ppval(s,t_fine),'r')

figure
hold on
plot(t_fine,ppval(p,t_fine),'k')
plot(t_fine,ppval(s,t_fine),'r')