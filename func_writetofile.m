function [io,OUTFILE] = func_writetofile(OUTFILE,io,OUT,grid,t,time,C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save model output to files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify variables to be written
if t==1
    io.varsout = [];
    io.unitsout = [];
    io.descout = [];
    io.typeout = [];
    if io.out_surface
        io.varsout{end+1} = 'smb';              io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Mass balance';
        io.varsout{end+1} = 'Tsurf';            io.unitsout{end+1} = 'K';               io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Surface temperature';
        io.varsout{end+1} = 'climT';            io.unitsout{end+1} = 'K';               io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Air temperature';
        io.varsout{end+1} = 'climP';            io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Precipitation';
        io.varsout{end+1} = 'climC';            io.unitsout{end+1} = 'fraction';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Cloud cover';
        io.varsout{end+1} = 'climRH';           io.unitsout{end+1} = 'fraction';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Relative humidity';
        io.varsout{end+1} = 'climWS';           io.unitsout{end+1} = 'm s-1';           io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Wind speed';
        io.varsout{end+1} = 'climPres';         io.unitsout{end+1} = 'Pa';              io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Air pressure';
        io.varsout{end+1} = 'climrain';         io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Rainfall';
        io.varsout{end+1} = 'climsnow';         io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Snowfall';
        io.varsout{end+1} = 'snowmass';         io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Snow mass';
        io.varsout{end+1} = 'mbal';             io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Cumulative mass balance';
        io.varsout{end+1} = 'melt';             io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Melt';
        io.varsout{end+1} = 'refr';             io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Refreezing';
        io.varsout{end+1} = 'runoff';           io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Runoff';
        io.varsout{end+1} = 'runoff_surf';      io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Surface runoff';
        io.varsout{end+1} = 'runoff_slush';     io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Slush runoff';
        io.varsout{end+1} = 'SWin';             io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Incoming SW radiation';
        io.varsout{end+1} = 'SWout';            io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Reflected SW radiation';
        io.varsout{end+1} = 'LWin';             io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Incoming LW radiation';
        io.varsout{end+1} = 'LWout';            io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Outgoing LW radiation';
        io.varsout{end+1} = 'SHF';              io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Sensible heat flux';
        io.varsout{end+1} = 'LHF';              io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Latent heat flux';
        io.varsout{end+1} = 'GHF';              io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Subsurface heat flux';
        io.varsout{end+1} = 'surfH';            io.unitsout{end+1} = 'm';               io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Surface height';
        io.varsout{end+1} = 'alb';              io.unitsout{end+1} = 'fraction';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Albedo';
    end
    if io.out_subsurface
        io.varsout{end+1} = 'subD';             io.unitsout{end+1} = 'kg m^{-3}';       io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Density';
        io.varsout{end+1} = 'subT';             io.unitsout{end+1} = 'K';               io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Temperature';         
        io.varsout{end+1} = 'subS';             io.unitsout{end+1} = 'mm w.e.';         io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Slush water content'; 
        io.varsout{end+1} = 'subW';             io.unitsout{end+1} = 'mm w.e.';         io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Irreducible water';   
        io.varsout{end+1} = 'subZ';             io.unitsout{end+1} = 'm';               io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Layer thickness';
    end
end


%% Update OUTFILE.TEMP with variables to be stored
fn = io.varsout;
for i=1:numel(fn)
    temp_long = eval(['OUT.' fn{i}]);
    if t==1
        OUTFILE.TEMP.(fn{i}) = zeros(size(temp_long));
    end
    if strcmp(io.typeout{i},'sample')
        if mod(t+floor(io.freqout/2),io.freqout)==0
            OUTFILE.TEMP.(fn{i}) = temp_long';
        end
    elseif strcmp(io.typeout{i},'mean')
        OUTFILE.TEMP.(fn{i}) = OUTFILE.TEMP.(fn{i}) + temp_long/io.freqout;
    elseif strcmp(io.typeout{i},'sum')
        OUTFILE.TEMP.(fn{i}) = OUTFILE.TEMP.(fn{i}) + temp_long;
    end
end

%% Save output to binary files
if t==1 
    if ~exist(io.outdir, 'dir')
       mkdir(io.outdir)
    end
    for i=1:length(fn)
        io.fid(i) = fopen([io.outdir 'OUT_' fn{i} '.bin'], 'w');
    end
end

if mod(t,io.freqout)==0
    for i=1:length(fn)
        OUTFILE.(fn{i}) = OUTFILE.TEMP.(fn{i});
        fwrite(io.fid(i),OUTFILE.(fn{i}),'real*4');
        OUTFILE.TEMP.(fn{i}) = 0.0;
    end
end

if t==time.tn
    for i=1:length(fn)
        fclose(io.fid(i));
    end
end

%% Save run info to file
if t==time.tn
    runinfo.grid = grid;
    runinfo.time = time;
    runinfo.IOout = io;
    runinfo.Cout = C;
    cd(io.outdir);
    save(io.infofile,'-struct','runinfo');
    cd(io.homedir);
end

end