%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNOW AND FIRN MODEL:
%%% - snow fall & riming
%%% - melting & sublimation
%%% - gravitational densification
%%% - heat diffusion
%%% - liquid water percolation, refreezing & storage
%%% - slush runoff and storage
%%% - refreezing stored slush water
%%% - refreezing stored irreducible water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT] = TIME_snow_model(C,OUT,IN,dt,grid,phys)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNOWFALL AND DEPOSITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fresh snow density [SOURCE: Kampenhout et al. (2017)]
if strcmp(phys.snow_compaction,'firn+snow')
    OUT.Dfreshsnow_T(IN.T>C.T0+2) = 50+1.7*17^(3/2);
    OUT.Dfreshsnow_T(IN.T<=C.T0+2 & IN.T>C.T0-15) = 50+1.7.*(IN.T(...
        IN.T<=C.T0+2 & IN.T>C.T0-15)-C.T0+15).^(3/2);
    OUT.Dfreshsnow_T(IN.T<=C.T0-15) = -3.8328*(IN.T(IN.T<=C.T0-15)-C.T0)-...
        0.0333*(IN.T(IN.T<=C.T0-15)-C.T0).^2;
    OUT.Dfreshsnow_W = 266.86.*(0.5.*(1+tanh(IN.WS./5))).^8.8;
    OUT.Dfreshsnow = OUT.Dfreshsnow_T(:) + OUT.Dfreshsnow_W(:);
else
    OUT.Dfreshsnow(:) = C.Dfreshsnow;
end

% Update layer depths & properties
shift_snowfall = IN.snow.*C.Dwater./OUT.Dfreshsnow;
shift_riming = OUT.moist_deposition.*C.Dwater./OUT.Dfreshsnow;
shift_tot(:,1) = shift_snowfall + shift_riming;
OUT.surfH(:,1) = OUT.surfH(:,1) + shift_tot(:,1);
runoff_irr_deep = zeros(grid.gpsum,1);
nl = grid.nl;

while any(shift_tot>0)
    
    shift = min(shift_tot,grid.max_subZ);
    shift_tot = shift_tot - shift;
        
    subT_old = OUT.subT;
    subD_old = OUT.subD;
    subW_old = OUT.subW;
    subZ_old = OUT.subZ;
    
    i_noshift = find(subZ_old(:,1)+shift<=grid.max_subZ);
    i_shift = find(subZ_old(:,1)+shift>grid.max_subZ);

    OUT.subZ(i_noshift,1) = subZ_old(i_noshift,1) + shift(i_noshift,1);
    OUT.subT(i_noshift,1) = subT_old(i_noshift,1).*subZ_old(i_noshift,1)...
        ./OUT.subZ(i_noshift,1) + ...
        OUT.Tsurf(i_noshift,1).*shift(i_noshift,1)./OUT.subZ(i_noshift,1);
    OUT.subD(i_noshift,1) = subD_old(i_noshift,1).*subZ_old(i_noshift,1)...
        ./OUT.subZ(i_noshift,1) + ...
        OUT.Dfreshsnow(i_noshift,1).*shift(i_noshift,1)./OUT.subZ(...
        i_noshift,1);
    OUT.subW(i_noshift,1) = subW_old(i_noshift,1);

    OUT.subZ(i_shift,3:nl) = subZ_old(i_shift,2:nl-1);
    OUT.subT(i_shift,3:nl) = subT_old(i_shift,2:nl-1);
    OUT.subD(i_shift,3:nl) = subD_old(i_shift,2:nl-1);
    OUT.subW(i_shift,3:nl) = subW_old(i_shift,2:nl-1);
    OUT.subZ(i_shift,2) = grid.max_subZ;
    OUT.subZ(i_shift,1) = (subZ_old(i_shift,1)+shift(i_shift,1)) - ...
        grid.max_subZ;
    OUT.subT(i_shift,2) = subT_old(i_shift,1).*subZ_old(i_shift,1)./...
        OUT.subZ(i_shift,2) + ...
        OUT.Tsurf(i_shift,1).*(OUT.subZ(i_shift,2)-subZ_old(i_shift,1))...
        ./OUT.subZ(i_shift,2);
    OUT.subT(i_shift,1) = OUT.Tsurf(i_shift,1);
    OUT.subD(i_shift,2) = subD_old(i_shift,1).*subZ_old(i_shift,1)...
        ./OUT.subZ(i_shift,2) + ...
        OUT.Dfreshsnow(i_shift,1).*(OUT.subZ(i_shift,2)-subZ_old(...
        i_shift,1))./OUT.subZ(i_shift,2);
    OUT.subD(i_shift,1) = OUT.Dfreshsnow(i_shift,1);
    OUT.subW(i_shift,2) = subW_old(i_shift,1);
    OUT.subW(i_shift,1) = 0.0;
    runoff_irr_deep(i_shift,1) = runoff_irr_deep(i_shift,1) + ...
        subW_old(i_shift,nl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MELT AND SUBLIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update layer depths & properties
OUT.sumWinit = sum(OUT.subW,2);
mass_removed = (OUT.melt + OUT.moist_sublimation)*1d3;
mass_layer = OUT.subD.*OUT.subZ;
n = 0;
while any(mass_removed>0)
    n = n+1;
    cond1 = mass_removed(:,1)>mass_layer(:,n);
    cond2 = ~cond1 & mass_removed(:,1)>0;
    mass_removed(cond1,1) = mass_removed(cond1,1) - OUT.subD(cond1,n).*...
        OUT.subZ(cond1,n);
    shift_tot(cond1,1) = shift_tot(cond1,1) - OUT.subZ(cond1,n);
    shift_tot(cond2,1) = shift_tot(cond2,1) - mass_removed(cond2,1)./...
        mass_layer(cond2,n).*OUT.subZ(cond2,n);
    mass_removed(cond2,1) = 0.0;
end

while any(shift_tot<0)
    shift = max(shift_tot,-OUT.subZ(:,2));
    shift_tot = shift_tot - shift;
    
    OUT.surfH(:,1) = OUT.surfH(:,1) + shift(:,1); 
    
    subT_old = OUT.subT;
    subD_old = OUT.subD;
    subW_old = OUT.subW;
    subZ_old = OUT.subZ;
    
    i_noshift = find(subZ_old(:,1)+shift>1d-17);
    i_shift = find(subZ_old(:,1)+shift<=1d-17);
    
    OUT.subZ(i_noshift,1) = subZ_old(i_noshift,1) + shift(i_noshift,1);
    OUT.subT(i_noshift,1) = subT_old(i_noshift,1);
    OUT.subD(i_noshift,1) = subD_old(i_noshift,1);
    temp = OUT.subZ(i_noshift,1)./subZ_old(i_noshift,1);
    OUT.subW(i_noshift,1) = subW_old(i_noshift,1).*temp;
    
    OUT.subZ(i_shift,2:nl-1) = subZ_old(i_shift,3:nl);
    OUT.subT(i_shift,2:nl-1) = subT_old(i_shift,3:nl);
    OUT.subD(i_shift,2:nl-1) = subD_old(i_shift,3:nl);
    OUT.subW(i_shift,2:nl-1) = subW_old(i_shift,3:nl);
    OUT.subZ(i_shift,1) = subZ_old(i_shift,1) + subZ_old(i_shift,2) + ...
        shift(i_shift,1);
    OUT.subT(i_shift,1) = subT_old(i_shift,2);
    OUT.subD(i_shift,1) = subD_old(i_shift,2);
    temp = OUT.subZ(i_shift,1)./subZ_old(i_shift,2);
    OUT.subW(i_shift,1) = subW_old(i_shift,2).*temp;
    for n=1:length(i_shift)
        if grid.doubledepth==1
            OUT.subZ(i_shift(n),nl) = 2.0^length(grid.split)*grid.max_subZ;
        else
            OUT.subZ(i_shift(n),nl) = grid.max_subZ;
        end
        OUT.subD(i_shift(n),nl) = subD_old(i_shift(n),nl);
    end
    OUT.subT(i_shift,nl) =  subT_old(i_shift,nl);
    OUT.subW(i_shift,nl) = 0.0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNOW COMPACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subD_old = OUT.subD;
subZ_old = OUT.subZ;
mliqmax = zeros(grid.gpsum,nl);

% Densification of firn [SOURCE: Ligtenberg et al. (2011)]
if strcmp(phys.snow_compaction,'firn_only')
    cond = true(grid.gpsum,nl);
elseif strcmp(phys.snow_compaction,'firn+snow')
    cond = OUT.subD>=C.Dfirn;
else
    error('phys.snow_compaction not set correctly!');
end

OUT.subTmean = OUT.subTmean .* (1-dt./C.yeardays) + OUT.subT .* dt./...
    C.yeardays;

grav_const = zeros(grid.gpsum,nl);
grav_const(cond) =  0.07 .* max(1.435-0.151*IN.logyearsnow(cond),0.25) ...
    .* (OUT.subD(cond)<550)+ 0.03 .* max(2.366-0.293*IN.logyearsnow(...
    cond),0.25) .* (OUT.subD(cond)>=550);   

OUT.subD(cond) = OUT.subD(cond) + dt./C.yeardays.*grav_const(cond).*...
    IN.yearsnow(cond).*C.g.*(C.Dice-OUT.subD(cond)).*exp(-C.Ec./(C.rd.*...
    OUT.subT(cond)) + C.Eg./(C.rd.*OUT.subTmean(cond)));
     
OUT.Dens_firn = zeros(grid.gpsum,nl);     
OUT.Dens_firn(cond) = dt./C.yeardays.*grav_const(cond).*IN.yearsnow(...
    cond).*C.g.*(C.Dice-OUT.subD(cond)).*exp(-C.Ec./(C.rd.*OUT.subT(...
    cond)) + C.Eg./(C.rd.*OUT.subTmean(cond)));
     
% Densification of seasonal snow [SOURCE: Kampenhout et al. (2017)]
if strcmp(phys.snow_compaction,'firn_only')
    cond = false(grid.gpsum,nl);
elseif strcmp(phys.snow_compaction,'firn+snow')
    cond = OUT.subD<C.Dfirn;
else
    error('phys.snow_compaction not set correctly!');
end

CC_tap = 175;
CC1 = zeros(grid.gpsum,nl);
CC2 = zeros(grid.gpsum,nl);
CC1(OUT.subD<CC_tap) = 1;
CC1(OUT.subD>=CC_tap) = exp(-0.046.*(OUT.subD(OUT.subD>=CC_tap)-CC_tap));
CC2(OUT.subW==0) = 1;
CC2(OUT.subW~=0) = 2;
CC3 = 2.777d-6;
CC4 = 0.04;

OUT.subD(cond) = OUT.subD(cond) + dt*C.dayseconds.*OUT.subD(cond).*CC3.* ...
    CC2(cond).*CC1(cond).*exp(-CC4.*(C.T0-OUT.subT(cond)));
OUT.subD(cond) = min(OUT.subD(cond),C.Dice);

OUT.Dens_destr_metam = zeros(grid.gpsum,nl);   
OUT.Dens_destr_metam(cond) = dt*C.dayseconds.*OUT.subD(cond).*CC3.*CC2(cond)...
    .*CC1(cond).*exp(-CC4.*(C.T0-OUT.subT(cond)));

CC5 = 0.1;
CC6 = 0.023;
CC7 = 4.0.*7.62237d6.*OUT.subD./358.0.*1./(1+60.*OUT.subW./(C.Dwater.*...
    OUT.subZ));
Psload = zeros(grid.gpsum,nl);
Visc = zeros(grid.gpsum,nl);
temp = cumsum(OUT.subD.*OUT.subZ,2) - 0.5*OUT.subD.*OUT.subZ;
Psload(cond) = temp(cond); 
Visc(cond) = CC7(cond).*exp(CC5.*(C.T0-OUT.subT(cond))+CC6.*...
    OUT.subD(cond));

OUT.subD(cond) = OUT.subD(cond) + dt*C.dayseconds.*OUT.subD(cond).*...
    Psload(cond)./Visc(cond);
OUT.subD(cond) = min(OUT.subD(cond),C.Dice);

OUT.Dens_overb_pres = zeros(grid.gpsum,nl); 
OUT.Dens_overb_pres(cond) = dt*C.dayseconds.*OUT.subD(cond).*Psload(cond)...
    ./Visc(cond);

MO = -0.069+0.66.*(1.25-0.0042.*(max(OUT.subD,50)-50));
SI = -2.868.*exp(-0.085.*IN.WS)+1+MO;
cond_SD = SI>0;
z_i = zeros(grid.gpsum,nl);
z_i(:,2:nl) = cumsum(OUT.subZ(:,1:nl-1).*(3.25-SI(:,1:nl-1)),2);
gamma_drift = max(0,SI.*exp(-z_i./0.1));
tau = 48*2*3600;
tau_i = tau./gamma_drift;

OUT.subD(cond & cond_SD) = OUT.subD(cond & cond_SD) + dt*3600*24.*max(...
    350-OUT.subD(cond & cond_SD),0)./tau_i(cond & cond_SD);
OUT.subD(cond & cond_SD) = min(OUT.subD(cond & cond_SD),C.Dice);

OUT.Dens_drift = zeros(grid.gpsum,nl); 
OUT.Dens_drift(cond & cond_SD) = dt*C.dayseconds.*max(350-OUT.subD(cond & ...
    cond_SD),0)./tau_i(cond & cond_SD);

% Update layer depths & properties
cond = OUT.subD<C.Dice;

OUT.subZ(cond) = subZ_old(cond).*subD_old(cond)./OUT.subD(cond);

mliqmax(cond) =  OUT.subD(cond).*OUT.subZ(cond).*0.0143.*exp(3.3.*(...      % SOURCE: Schneider and Jansson (2004)
    C.Dice-OUT.subD(cond))./C.Dice)./(1-0.0143.*exp(3.3.*(C.Dice-...
    OUT.subD(cond))./C.Dice))...
    .*0.05.*min(C.Dice-OUT.subD(cond),20);
OUT.subW = min(mliqmax,OUT.subW);

shift(:,1) = sum(OUT.subZ,2)-sum(subZ_old,2);
OUT.surfH(:,1) = OUT.surfH(:,1) + shift(:,1);

runoff_irr = OUT.sumWinit - sum(OUT.subW,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HEAT CONDUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dz1 = ((OUT.subZ(:,1)+0.5.*OUT.subZ(:,2)).^2);
dz2 = (0.5.*(OUT.subZ(:,3:nl)+OUT.subZ(:,2:nl-1)).^2);
kk = zeros(grid.gpsum,nl);
c_eff = zeros(grid.gpsum,nl);
kdTdz = zeros(grid.gpsum,nl);

% effective conductivity
kk = 0.138-1.01d-3.*OUT.subD+3.233d-6.*OUT.subD.^2;                         % SOURCE: Sturm et al. (1997)
 
% specific heat capacity
c_eff = OUT.subD .* (152.2+7.122.*OUT.subT);                                % SOURCE: Yen (1981)                              
    
% adaptive time-stepping to assure stability (CFL condition)
z_temp = OUT.subZ(:,2:end);
c_eff_temp = c_eff(:,1:end);
kk_temp = kk(:,1:end);
dt_stab = 0.5*min(c_eff_temp,[],2).*min(z_temp,[],2).^2./max(...
    kk_temp,[],2)/C.dayseconds;

% update temperatures after heat conduction / diffusion
tt = zeros(grid.gpsum,1);
while any(tt<dt)
    subT_old = OUT.subT;
    
    dt_temp = min(dt_stab,-tt+dt);

    tt = tt + dt_temp;
    
    cond_dt = dt_temp>0;
    
    kdTdz(cond_dt,2) = (kk(cond_dt,1).*OUT.subZ(cond_dt,1)+0.5.*kk(...      % subsurface heat flux (at upper boundary of layer 2)
        cond_dt,2).*OUT.subZ(cond_dt,2)).*(subT_old(cond_dt,2)-...
        OUT.Tsurf(cond_dt,1)) ./ dz1(cond_dt,1);
    kdTdz(cond_dt,3:nl) = (kk(cond_dt,2:nl-1).*OUT.subZ(cond_dt,2:nl-1)+... % subsurface heat flux (at upper boundary of layer 3 ... nl) 
        kk(cond_dt,3:nl).*OUT.subZ(cond_dt,3:nl)).*(subT_old(...
        cond_dt,3:nl)-subT_old(cond_dt,2:nl-1)) ./ dz2(cond_dt,:);

    OUT.subT(cond_dt,2) = subT_old(cond_dt,2) + C.dayseconds.*dt_temp(...   % temperature 2nd layer 
        cond_dt,1).*(kdTdz(cond_dt,3)-kdTdz(cond_dt,2))./(c_eff(...
        cond_dt,2).*(0.5*OUT.subZ(cond_dt,1)+0.5*OUT.subZ(cond_dt,2)+...
        0.25*OUT.subZ(cond_dt,3)));
    OUT.subT(cond_dt,3:nl-1) = subT_old(cond_dt,3:nl-1) + C.dayseconds.*... % temperature layer 3 ... nl-1
        dt_temp(cond_dt,1).*(kdTdz(cond_dt,4:nl)-kdTdz(cond_dt,3:nl-1))...
        ./(c_eff(cond_dt,3:nl-1).*(0.25*OUT.subZ(cond_dt,2:nl-2)+0.5*...
        OUT.subZ(cond_dt,3:nl-1)+0.25*OUT.subZ(cond_dt,4:nl)));
                    
    OUT.subT(cond_dt,nl) = subT_old(cond_dt,nl) + C.dayseconds.*dt_temp(... % temperature bottom layer
        cond_dt,1).*(C.geothermal_flux-kdTdz(cond_dt,nl))./(c_eff(...
        cond_dt,nl).*(0.25*OUT.subZ(cond_dt,nl-1)+0.75*OUT.subZ(...
        cond_dt,nl)));
end
OUT.subT(:,1) = OUT.Tsurf(:,1) + (OUT.subT(:,2)-OUT.Tsurf(:,1))./...        % temperature first layer
    (OUT.subZ(:,1)+0.5*OUT.subZ(:,2)).*0.5.*OUT.subZ(:,1);
OUT.subT(OUT.subT>C.T0) = C.T0;
OUT.subCeff = c_eff;
OUT.subK = kk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WATER PERCOLATION, REFREEZING AND IRREDUCIBLE WATER STORAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subW_old = OUT.subW;

% Water input
avail_W =   OUT.melt*1d3 + ...                                              % melt water
            IN.rain.*1d3 + ...                                              % rain fall         
            (OUT.moist_condensation-OUT.moist_evaporation).*1d3;            % condensation or evaporation                                   
avail_W = max(avail_W,0);

% Check conditions for refreezing and irreducible water storage
OUT.cpi = 152.2+7.122.*OUT.subT;
c1 = OUT.cpi.*OUT.subD.*OUT.subZ.*(C.T0-OUT.subT)./C.Lm;                    % cold content limit on refreezing
c2 = OUT.subZ.*(1-OUT.subD./C.Dice).*C.Dice;                                % density limit on refreezing
cond1 = c1 >= c2;
cond2 = ~cond1;
Wlim = zeros(grid.gpsum,nl);                                                % maximum potential for refreezing
Wlim(cond1) = c2(cond1);
Wlim(cond2) = c1(cond2);
Wlim = max(Wlim,0);
mliqmax = zeros(grid.gpsum,nl);                                             % maximum irreducible water storage
noice = OUT.subD<C.Dice-1;
mliqmax(noice) =  OUT.subD(noice).*OUT.subZ(noice).*0.0143.*exp(3.3.*...
    (C.Dice-OUT.subD(noice))./C.Dice)./(1-0.0143.*exp(3.3.*...
    (C.Dice-OUT.subD(noice))./C.Dice)).*0.05.*min(C.Dice-...
    OUT.subD(noice),20);
Wirr = mliqmax-subW_old;
RP = zeros(grid.gpsum,nl);                                                 
leftW = zeros(grid.gpsum,1);
avail_W_loc = zeros(grid.gpsum,1);

% Water percolation
z0 = C.perc_depth;
zz = cumsum(OUT.subZ,2)-0.5*OUT.subZ;
carrot = zeros(grid.gpsum,nl);      
if  strcmp(phys.percolation,'bucket')                                       % all water added at surface
    carrot(:,1) = 1;
elseif   strcmp(phys.percolation,'normal')                                  % normal percolation distribution
    carrot = 2*exp( - zz.^2/2/(z0/3)^2 )/(z0/3)/sqrt(2*pi);
elseif   strcmp(phys.percolation,'linear')                                  % linear percolation distribution
    carrot = 2*(z0 - zz) / z0.^2;
    carrot = max(carrot,0);                                               
elseif   strcmp(phys.percolation,'uniform')                                 % uniform percolation distribution
    [~, ind] = min(abs(zz - z0));
    carrot(1:ind) = 1/z0; clear ind
else
    error('phys.percolation is not set correctly!');
end
carrot = carrot.*OUT.subZ;
carrot = carrot./sum(carrot,2);
carrot = carrot.*avail_W;

% Refreezing and irreducible water storage per layer
for n=1:nl
    avail_W_loc = avail_W_loc + carrot(:,n);                                % available water per layer
    
    cond1 = avail_W_loc>Wlim(:,n);
    
    RP(cond1,n) = Wlim(cond1,n);                                            % more water than refreezing limit
    leftW(cond1,1) = (avail_W_loc(cond1,1)-Wlim(cond1,n));
    OUT.subW(cond1,n) = subW_old(cond1,n) + min(leftW(cond1,1),...
        Wirr(cond1,n));
    
    RP(~cond1,n) = avail_W_loc(~cond1);                                     % less water than refreezing limit
    OUT.subW(~cond1,n) = subW_old(~cond1,n);
    
    avail_W_loc = avail_W_loc - RP(:,n) - (OUT.subW(:,n)-subW_old(:,n));
    
    OUT.subT(:,n) = OUT.subT(:,n) + C.Lm.*RP(:,n)./(OUT.subD(:,n).*...      % update temperature after refreezing
        OUT.cpi(:,n).*OUT.subZ(:,n));
    OUT.subD(:,n) = OUT.subD(:,n) + RP(:,n)./OUT.subZ(:,n);                 % update density after refreezing
end
avail_W = avail_W_loc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SLUSH WATER STORAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Available pore space for storing slush water
slushspace = max(OUT.subZ.*(1-OUT.subD./C.Dice).*C.Dwater-OUT.subW,0);      % available pore volume to store slush wate
total_slushspace = sum(slushspace,2);

% Calculate surface runoff
avail_W = avail_W + sum(OUT.subS,2);
runoff_surface = max(avail_W - total_slushspace,0);                         % runoff of excess water at surface

% Update slush water content after new water input and runoff
avail_S = min(avail_W,total_slushspace);
runoff_slush = avail_S - 1.0/(1.0+dt/C.Trunoff).*avail_S;                   % runoff of slush water
avail_S = 1.0/(1.0+dt/C.Trunoff).*avail_S;
avail_S(avail_S<1d-25) = 0.0;
OUT.subS = zeros(grid.gpsum,nl);
for n=nl:-1:1
    cond1 = avail_S>slushspace(:,n);
    OUT.subS(cond1,n) = slushspace(cond1,n);
    OUT.subS(~cond1,n) = avail_S(~cond1);
    avail_S = avail_S - OUT.subS(:,n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REFREEZING OF STORED WATER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refreezing of slush water
OUT.cpi = 152.2+7.122.*OUT.subT;                                            % specific heat capacity
c1 = OUT.cpi.*OUT.subD.*OUT.subZ.*(C.T0-OUT.subT)./C.Lm;                    % cold content limit on refreezing
c2 = OUT.subZ.*(1-OUT.subD./C.Dice).*C.Dice;                                % density limit on refreezing
cond1 = c1 >= c2;
cond2 = ~cond1;
Wlim = zeros(grid.gpsum,nl);                                                % maximum potential for refreezing
Wlim(cond1) = c2(cond1);
Wlim(cond2) = c1(cond2);
RS = zeros(grid.gpsum,nl);                                                  % amount of slush water refreezing

for n=nl:-1:1
    if (any(OUT.subS(:,n)>0) && any(OUT.subT(:,n)<C.T0))
        avail_W = OUT.subS(:,n);
        Wlim_loc = Wlim(:,n);

        cond = avail_W>Wlim_loc;
        RS(cond,n) = Wlim_loc(cond);
        RS(~cond,n) = avail_W(~cond);

        OUT.subS(:,n) = OUT.subS(:,n) - RS(:,n);                            % update slush water content after refreezing
        OUT.subT(:,n) = OUT.subT(:,n) + C.Lm.*RS(:,n)./(OUT.subD(:,n)...    % update temperature after refreezing
            .*OUT.cpi(:,n).*OUT.subZ(:,n));
        OUT.subD(:,n) = OUT.subD(:,n) + RS(:,n)./OUT.subZ(:,n);             % update density after refreezing
    end
end

% Refreezing of irreducible water
OUT.cpi = 152.2+7.122.*OUT.subT;                                            % specific heat capacity
c1 = OUT.cpi.*OUT.subD.*OUT.subZ.*(C.T0-OUT.subT)./C.Lm;                    % cold content limit on refreezing
c2 = OUT.subZ.*(1-OUT.subD./C.Dice).*C.Dice;                                % density limit on refreezing
cond1 = c1 >= c2;
cond2 = ~cond1;
Wlim = zeros(grid.gpsum,nl);                                                % maximum potential for refreezing                        
Wlim(cond1) = c2(cond1);
Wlim(cond2) = c1(cond2);
RI = zeros(grid.gpsum,nl);                                                  % amount of irreducible water refreezing

for n=nl:-1:1
    if (any(OUT.subW(:,n)>0) && any(OUT.subT(:,n)<C.T0))
        avail_W = OUT.subW(:,n);                                            % refreeze irreducible water
        cond1 = avail_W>Wlim(:,n);
        RI(cond1,n) = Wlim(cond1,n);
        RI(~cond1,n) = avail_W(~cond1);

        OUT.subW(:,n) = OUT.subW(:,n) - RI(:,n);                            % update water content after refreezing
        OUT.subT(:,n) = OUT.subT(:,n) + C.Lm.*RI(:,n)./(OUT.subD(:,n)...    % update temperature after refreezing
            .*OUT.cpi(:,n).*OUT.subZ(:,n));
        OUT.subD(:,n) = OUT.subD(:,n) + RI(:,n)./OUT.subZ(:,n);             % update density after refreezing
    end
end

OUT.refr = 1d-3.*(sum(RP,2) +  sum(RS,2) +  sum(RI,2));                     % total refreezing
OUT.refr_P = 1d-3.*sum(RP,2);                                               % refreezing of percolating water
OUT.refr_S = 1d-3.*sum(RS,2);                                               % refreezing of slush water
OUT.refr_I = 1d-3.*sum(RI,2);                                               % refreezing of irreducible water
OUT.slushw = sum(OUT.subS,2);                                               % total stored slush water
OUT.irrw = sum(OUT.subW,2);                                                 % total stored irreducible water

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VERTICAL LAYER MERGING AND SPLITTING (IN CASE OF LAYER DOUBLING)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if grid.doubledepth
    for n=1:length(grid.split)
        % Merge layers (accumulation case)
        split = grid.split(n);
        cond = OUT.subZ(:,split)<=(2.0^(n-1))*grid.max_subZ & ...
            grid.mask==1;
        subT_old = OUT.subT;
        subD_old = OUT.subD;
        subW_old = OUT.subW;
        subZ_old = OUT.subZ;
        subS_old = OUT.subS;
        OUT.subZ(cond,split-1) = subZ_old(cond,split-1) + ...
            subZ_old(cond,split);
        OUT.subW(cond,split-1) = subW_old(cond,split-1) + ...
            subW_old(cond,split);
        OUT.subS(cond,split-1) = subS_old(cond,split-1) + ...
            subS_old(cond,split);
        OUT.subD(cond,split-1) = (subZ_old(cond,split-1).*...
            subD_old(cond,split-1) + subZ_old(cond,split).*...
            subD_old(cond,split))./(subZ_old(cond,split-1) + ...
            subZ_old(cond,split));
        OUT.subT(cond,split-1) = (subZ_old(cond,split-1).*...
            subT_old(cond,split-1) + subZ_old(cond,split).*...
            subT_old(cond,split))./(subZ_old(cond,split-1) + ...
            subZ_old(cond,split));
        OUT.subZ(cond,split:nl-1) = subZ_old(cond,split+1:nl);
        OUT.subW(cond,split:nl-1) = subW_old(cond,split+1:nl);
        OUT.subS(cond,split:nl-1) = subS_old(cond,split+1:nl);
        OUT.subD(cond,split:nl-1) = subD_old(cond,split+1:nl);
        OUT.subT(cond,split:nl-1) = subT_old(cond,split+1:nl);
        OUT.subZ(cond,nl) = 2.0^length(grid.split)*grid.max_subZ;
        OUT.subT(cond,nl) = 2.0*subT_old(cond,nl) - subT_old(cond,nl-1);
        OUT.subD(cond,nl) = subD_old(cond,nl);
        OUT.subW(cond,nl) = 0.0;
        OUT.subS(cond,nl) = 0.0;
        
        % Split layers (ablation case)
        cond = OUT.subZ(:,split-2)>(2.0^(n-1))*grid.max_subZ & ...
            grid.mask==1;
        subT_old = OUT.subT;
        subD_old = OUT.subD;
        subW_old = OUT.subW;
        subZ_old = OUT.subZ;
        subS_old = OUT.subS;
        OUT.subZ(cond,split-2) = 0.5*subZ_old(cond,split-2);
        OUT.subW(cond,split-2) = 0.5*subW_old(cond,split-2);
        OUT.subS(cond,split-2) = 0.5*subS_old(cond,split-2);
        OUT.subT(cond,split-2) = subT_old(cond,split-2);
        OUT.subD(cond,split-2) = subD_old(cond,split-2);
        OUT.subZ(cond,split-1) = 0.5*subZ_old(cond,split-2);
        OUT.subW(cond,split-1) = 0.5*subW_old(cond,split-2);
        OUT.subS(cond,split-1) = 0.5*subS_old(cond,split-2);
        OUT.subT(cond,split-1) = subT_old(cond,split-2);
        OUT.subD(cond,split-1) = subD_old(cond,split-2);
        OUT.subZ(cond,split:nl) = subZ_old(cond,split-1:nl-1);
        OUT.subW(cond,split:nl) = subW_old(cond,split-1:nl-1);
        OUT.subS(cond,split:nl) = subS_old(cond,split-1:nl-1);
        OUT.subT(cond,split:nl) = subT_old(cond,split-1:nl-1);
        OUT.subD(cond,split:nl) = subD_old(cond,split-1:nl-1);
        runoff_irr_deep(cond,1) = runoff_irr_deep(cond,1) + ...
            subW_old(cond,nl);
        runoff_slush(cond,1) = runoff_slush(cond,1) + subS_old(cond,nl);     
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUNOFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT.runoff_irr_deep_mean = OUT.runoff_irr_deep_mean .* ...                  % runoff of irreducible water below the model bottom
    (1-dt./C.yeardays) + runoff_irr_deep .* dt./C.yeardays;         
OUT.runoff = 1d-3.*(runoff_surface + runoff_slush + runoff_irr ...          % total runoff
    + OUT.runoff_irr_deep_mean);
OUT.runoff_surf = 1d-3.*runoff_surface;                                     % surface runoff [in m w.e. per timestep]
OUT.runoff_slush = 1d-3.*runoff_slush;                                      % slush runoff [in m w.e. per timestep]
OUT.runoff_irr = 1d-3.*runoff_irr;                                          % irreducible water runoff within domain [in m w.e. per timestep]
OUT.runoff_irr_deep = 1d-3.*OUT.runoff_irr_deep_mean;                       % irreducible water runoff below domain [in m w.e. per timestep]

end

