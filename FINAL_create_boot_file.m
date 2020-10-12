%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CREATE RESTART FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FINAL_create_boot_file(OUT,io)

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

if ~exist(io.rebootdir, 'dir')
    mkdir(io.rebootdir);
end
if (io.writebootfile)
    save([io.rebootdir io.bootfileout],'boot');
end

end







