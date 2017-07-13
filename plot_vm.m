m = 5;
l = 3;


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

%load in best-fit FRET sims
load(['FRET_interp_est_mean' num2str(well) '_final.mat'])

%fret as a function of time
FRET = @(t) ppval(p,t);

%points along m that we care about
n=5;
msamp = augknt([m0,m1,linspace(m0,m1,n)],2);
%create spline functions
v_spline = spmak(msamp,q_final{l}');
%evaluate V(m(t))
Vx_c = @(t) fnval(v_spline,FRET(t));

V_m = @(t) fnval(v_spline,t);

figure
% subplot(1,2,1)

plot(linspace(m0,m1,100),V_m(linspace(m0,m1,100)))

xlabel('m')
ylabel('v(m)')
title('Spline interpolation of v(m)')

exportfig(gcf,['FRET_speed_interp_' num2str(m) '.eps'],'color','rgb')

% subplot(1,2,2)
% 
% plot(tdata,Vx_c(tdata)./max(Vx_c(tdata)))
% hold on
% plot(tdata,FRET(tdata),'r')
