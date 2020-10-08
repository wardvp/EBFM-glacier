%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save output file for rebooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func_createbootfile(OUT,io)
boot = struct;
boot.subZ = OUT.subZ;
boot.subW = OUT.subW;
boot.subD = OUT.subD;
boot.subS = OUT.subS;
boot.subT = OUT.subT;
boot.subSS = OUT.subSS;
boot.subSOIL = OUT.subSOIL;
boot.subTmean = OUT.subTmean;
boot.mbal = OUT.mbal;
boot.mbal_stake = OUT.mbal_stake;
boot.snowmass = OUT.snowmass;
boot.Tsurf = OUT.Tsurf;
boot.ys = OUT.ys;
boot.timelastsnow = OUT.timelastsnow;
boot.alb_snow = OUT.alb_snow;

cd(io.rebootdir);
if (io.writebootfile)
    save(io.bootfileout,'boot');
end
cd(io.homedir);

end







