%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURES: DISTRIBUTED GRID
%%% Plots temporal mean spatial maps of a specified variable for a
%%% specified time-period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER SPECIFICATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdir = '..\Output\';                      % output directory
timedist_start = '1-Sep-2015 0:00';         % start date
timedist_end = '1-Sep-2016 0:00';           % start time
var = 'cmb_cumulative';                     % plot variable (choose from 
                                            %   list in io.varsout)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([outdir '\runinfo.mat']);

ind = find(strcmp(var,{IOout.varsout.varname}));

L = length(grid.mask);
time_start = time.ts;
time_end = time.te;
period = IOout.freqout*time.dt;
tvec = datenum(time_start):period:datenum(time_end);
T = (datenum(time_end) - datenum(time_start))/period;

varabr = IOout.varsout(ind).varname;
varname = IOout.varsout(ind).description;
units = IOout.varsout(ind).units;
type = IOout.varsout(ind).type;

taxis = datetime(datevec(datenum(timedist_start):IOout.freqout*time.dt:...
    datenum(timedist_end)+0.5*IOout.freqout*time.dt));
tind = [];
taxis_plot_datenum = [];
for tt=1:length(taxis)
    if min(abs(datenum(IOout.output_times)-datenum(taxis(tt))))<...
            IOout.freqout*time.dt
        [~,tind(end+1)] = min(abs(datenum(IOout.output_times)-...
            datenum(taxis(tt))));
        taxis_plot_datenum(end+1) = datenum(IOout.output_times(tind(end)));
    end
end
taxis_plot = datetime(datevec(taxis_plot_datenum));


%% Distributed plot (time-averaged)
A2D = nan(grid.Lx,grid.Ly);
A2D(grid.ind) = 0;

fid = fopen([outdir '/OUT_' varabr '.bin'],'rb');
count = 0;
for t=1:length(tind)
    fseek(fid,(tind(t)-1)*4*L,'bof');
    Atemp = fread(fid,L,'real*4','l');
    A2D(grid.ind) = A2D(grid.ind) + Atemp;
    count = count + 1;
    A2D(grid.mask_2D==0) = NaN;
end
fclose(fid);
if strcmp(type,'mean')
    A2D = A2D / count;
end

%% FIGURE
f = figure;
pcolor(grid.x_2D,grid.y_2D,A2D); shading flat;
xlabel('UTM Easting (m)');
ylabel('UTM Northing (m)');
c = colorbar;
ylabel(c,[varname ' (' units ')']);
title(varname);
axis equal;
xlim([min(grid.x) max(grid.x)]);
ylim([min(grid.y) max(grid.y)]);