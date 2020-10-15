%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WRITE MODEL OUTPUT TO FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [io,OUTFILE] = TIME_write_to_file(OUTFILE,io,OUT,grid,t,time,C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Specify variables to be written
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==1
    io.output_times = [];
    OUTFILE.varsout = cell(1,4);
    OUTFILE.varsout(1,:)     = {'cmb'           ;'m w.e.'   ;'sum'      ;'Climatic mass balance'};
    OUTFILE.varsout(end+1,:) = {'Tsurf'         ;'K'        ;'mean'     ;'Surface temperature'};
    OUTFILE.varsout(end+1,:) = {'climT'         ;'K'        ;'mean'     ;'Air temperature'};
    OUTFILE.varsout(end+1,:) = {'climP'         ;'m w.e.'   ;'sum'      ;'Precipitation'};
    OUTFILE.varsout(end+1,:) = {'climC'         ;'fraction' ;'mean'     ;'Cloud cover'};
    OUTFILE.varsout(end+1,:) = {'climRH'        ;'fraction' ;'mean'     ;'Relative humidity'};
    OUTFILE.varsout(end+1,:) = {'climWS'        ;'m s-1'    ;'mean'     ;'Wind speed'};
    OUTFILE.varsout(end+1,:) = {'climPres'      ;'Pa'       ;'mean'     ;'Air pressure'};
    OUTFILE.varsout(end+1,:) = {'climrain'      ;'m w.e.'   ;'sum'      ;'Rainfall'};
    OUTFILE.varsout(end+1,:) = {'climsnow'      ;'m w.e.'   ;'sum'      ;'Snowfall'};
    OUTFILE.varsout(end+1,:) = {'snowmass'      ;'m w.e.'   ;'mean'     ;'Snow mass'};
    OUTFILE.varsout(end+1,:) = {'cmb_cumulative';'m w.e.'   ;'mean'     ;'Cumulative mass balance'};
    OUTFILE.varsout(end+1,:) = {'melt'          ;'m w.e.'   ;'sum'      ;'Melt'};
    OUTFILE.varsout(end+1,:) = {'refr'          ;'m w.e.'   ;'sum'      ;'Refreezing'};
    OUTFILE.varsout(end+1,:) = {'runoff'        ;'m w.e.'   ;'sum'      ;'Runoff'};
    OUTFILE.varsout(end+1,:) = {'runoff_surf'   ;'m w.e.'   ;'sum'      ;'Surface runoff'};
    OUTFILE.varsout(end+1,:) = {'runoff_slush'  ;'m w.e.'   ;'sum'      ;'Slush runoff'};
    OUTFILE.varsout(end+1,:) = {'SWin'          ;'W m^{-2}' ;'mean'     ;'Incoming SW radiation'};
    OUTFILE.varsout(end+1,:) = {'SWout'         ;'W m^{-2}' ;'mean'     ;'Reflected SW radiation'};
    OUTFILE.varsout(end+1,:) = {'LWin'          ;'W m^{-2}' ;'mean'     ;'Incoming LW radiation'};
    OUTFILE.varsout(end+1,:) = {'LWout'         ;'W m^{-2}' ;'mean'     ;'Outgoing LW radiation'};
    OUTFILE.varsout(end+1,:) = {'SHF'           ;'W m^{-2}' ;'mean'     ;'Sensible heat flux'};
    OUTFILE.varsout(end+1,:) = {'LHF'           ;'W m^{-2}' ;'mean'     ;'Latent heat flux'};
    OUTFILE.varsout(end+1,:) = {'GHF'           ;'W m^{-2}' ;'mean'     ;'Subsurface heat flux'};
    OUTFILE.varsout(end+1,:) = {'surfH'         ;'m'        ;'sample'   ;'Surface height'};
    OUTFILE.varsout(end+1,:) = {'alb'           ;'fraction' ;'mean'     ;'Albedo'};
    OUTFILE.varsout(end+1,:) = {'shade'         ;'fraction' ;'mean'     ;'Shading (0=not shaded, 1=shaded)'};
    
    OUTFILE.varsout(end+1,:) = {'subD'          ;'kg m^{-3}';'sample'   ;'Density'};
    OUTFILE.varsout(end+1,:) = {'subT'          ;'K'        ;'sample'   ;'Temperature'};
    OUTFILE.varsout(end+1,:) = {'subS'          ;'mm w.e.'  ;'sample'   ;'Slush water content'};
    OUTFILE.varsout(end+1,:) = {'subW'          ;'mm w.e.'  ;'sample'   ;'Irreducible water'};
    OUTFILE.varsout(end+1,:) = {'subZ'          ;'m'        ;'sample'   ;'Layer thickness'};
    
    io.varsout = struct('varname',{OUTFILE.varsout{1:end,1}},...
        'units',{OUTFILE.varsout{1:end,2}},'type',{OUTFILE.varsout{1:end,3}},...
        'description',{OUTFILE.varsout{1:end,4}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update OUTFILE.TEMP with variables to be stored
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = OUTFILE.varsout;
for i=1:length(fn)
    temp_long = OUT.(fn{i,1});
    if mod(t-1,io.freqout)==0
        OUTFILE.TEMP.(fn{i,1}) = zeros(size(temp_long));
    end
    if strcmp(fn{i,3},'sample')
        if mod(t+floor(io.freqout/2),io.freqout)==0
            OUTFILE.TEMP.(fn{i,1}) = temp_long';
        end
    elseif strcmp(fn{i,3},'mean')
        OUTFILE.TEMP.(fn{i,1}) = OUTFILE.TEMP.(fn{i,1}) + temp_long/io.freqout;
    elseif strcmp(fn{i,3},'sum')
        OUTFILE.TEMP.(fn{i,1}) = OUTFILE.TEMP.(fn{i,1}) + temp_long;
    end
end
if mod(t-1,io.freqout)==0
    OUTFILE.output_time = time.TCUR/io.freqout;
else
    OUTFILE.output_time = OUTFILE.output_time + time.TCUR/io.freqout;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save output to binary files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==1 
    if ~exist(io.outdir, 'dir')
       mkdir(io.outdir)
    end
    for i=1:length(fn)
        io.fid(i) = fopen([io.outdir 'OUT_' fn{i,1} '.bin'], 'w');
    end
end

if mod(t,io.freqout)==0
    for i=1:length(fn)
        OUTFILE.(fn{i,1}) = OUTFILE.TEMP.(fn{i,1});
        fwrite(io.fid(i),OUTFILE.(fn{i,1}),'real*4');
    end
    io.output_times(end+1,1:6) = datevec(OUTFILE.output_time);
end

if t==time.tn
    for i=1:length(fn)
        fclose(io.fid(i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save run info to runinfo.mat in the output directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==time.tn
    io.output_times = datetime(io.output_times);
    runinfo.grid = grid;
    runinfo.time = time;
    runinfo.IOout = io;
    runinfo.Cout = C;
    save([io.outdir 'runinfo.mat'],'-struct','runinfo');
end

end