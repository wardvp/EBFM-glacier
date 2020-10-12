%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET TIME PARAMETERS & PRINT MODEL TIME TO SCREEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time] = TIME_print_time(t,time)

time.TS = datenum(time.ts);
time.TE = datenum(time.te);

time.TCUR = time.TS+(t-1)*time.dt;
time.TCUR_DT = datetime(datevec(time.TCUR));
time.TPREV = time.TS+(t-2)*time.dt;
time.TPREV_DT = datetime(datevec(time.TPREV));

tempdate = datevec(time.TCUR);

fprintf('%6s','Year: ');    fprintf('%3i',tempdate(1));
fprintf('%9s','   Month: ');  fprintf('%3i',tempdate(2));
fprintf('%7s','   Day: ');  fprintf('%3i',tempdate(3));
fprintf('%8s','   Hour: ');  fprintf('%3i\n',tempdate(4));

end



