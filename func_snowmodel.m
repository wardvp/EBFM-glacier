%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update density, temperature, water content of vertical model after:
%%% - snow fall & riming
%%% - melting & sublimation
%%% - gravitational densification
%%% - heat diffusion
%%% - liquid water percolation, refreezing & storage
%%% - slush runoff and storage
%%% - refreezing stored slush water
%%% - refreezing stored irreducible water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT] = func_snowmodel(C,OUT,IN,dt,grid,time,phys)

nl = grid.nl;
OUT.sumSinit = sum(OUT.subS,2);
OUT.subTmean = OUT.subTmean .* (1-dt./C.yeardays) + OUT.subT .* dt./C.yeardays;

%% Fresh snow density 
% Temperature-dependent contribution [Kampenhout et al. 2017]
OUT.Dfreshsnow_T(IN.T>C.T0+2) = 50+1.7*17^(3/2);
OUT.Dfreshsnow_T(IN.T<=C.T0+2 & IN.T>C.T0-15) = 50+1.7.*(IN.T(IN.T<=C.T0+2 & IN.T>C.T0-15)-C.T0+15).^(3/2);
OUT.Dfreshsnow_T(IN.T<=C.T0-15) = -3.8328*(IN.T(IN.T<=C.T0-15)-C.T0)-0.0333*(IN.T(IN.T<=C.T0-15)-C.T0).^2;

% Wind-dependent contribution [Kampenhout et al. 2017]
OUT.Dfreshsnow_W = 266.86.*(0.5.*(1+tanh(IN.WS./5))).^8.8;

OUT.Dfreshsnow = OUT.Dfreshsnow_T(:) + OUT.Dfreshsnow_W(:);

%% Update vertical layers after snow fall & riming
shift_snowfall = IN.snow.*C.Dwater./OUT.Dfreshsnow;
shift_riming = OUT.moist_deposition.*C.Dwater./OUT.Dfreshsnow;
shift_tot(:,1) = shift_snowfall + shift_riming;
OUT.surfH(:,1) = OUT.surfH(:,1) + shift_tot(:,1);
runoff_irr_deep = zeros(grid.gpsum,1);

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
    OUT.subT(i_noshift,1) = subT_old(i_noshift,1).*subZ_old(i_noshift,1)./OUT.subZ(i_noshift,1) + ...
        OUT.Tsurf(i_noshift,1).*shift(i_noshift,1)./OUT.subZ(i_noshift,1);
    OUT.subD(i_noshift,1) = subD_old(i_noshift,1).*subZ_old(i_noshift,1)./OUT.subZ(i_noshift,1) + ...
        OUT.Dfreshsnow(i_noshift,1).*shift(i_noshift,1)./OUT.subZ(i_noshift,1);
    OUT.subW(i_noshift,1) = subW_old(i_noshift,1);

    % shift layers down when first layer exceeds maximum thickness
    OUT.subZ(i_shift,3:nl) = subZ_old(i_shift,2:nl-1);
    OUT.subT(i_shift,3:nl) = subT_old(i_shift,2:nl-1);
    OUT.subD(i_shift,3:nl) = subD_old(i_shift,2:nl-1);
    OUT.subW(i_shift,3:nl) = subW_old(i_shift,2:nl-1);
    OUT.subZ(i_shift,2) = grid.max_subZ;
    OUT.subZ(i_shift,1) = (subZ_old(i_shift,1)+shift(i_shift,1)) - grid.max_subZ;
    OUT.subT(i_shift,2) = subT_old(i_shift,1).*subZ_old(i_shift,1)./OUT.subZ(i_shift,2) + ...
        OUT.Tsurf(i_shift,1).*(OUT.subZ(i_shift,2)-subZ_old(i_shift,1))./OUT.subZ(i_shift,2);
    OUT.subT(i_shift,1) = OUT.Tsurf(i_shift,1);
    OUT.subD(i_shift,2) = subD_old(i_shift,1).*subZ_old(i_shift,1)./OUT.subZ(i_shift,2) + ...
        OUT.Dfreshsnow(i_shift,1).*(OUT.subZ(i_shift,2)-subZ_old(i_shift,1))./OUT.subZ(i_shift,2);
    OUT.subD(i_shift,1) = OUT.Dfreshsnow(i_shift,1);
    OUT.subW(i_shift,2) = subW_old(i_shift,1);
    OUT.subW(i_shift,1) = 0.0;
    runoff_irr_deep(i_shift,1) = runoff_irr_deep(i_shift,1) + subW_old(i_shift,nl);
end

%% Update vertical layers after melt & sublimation
OUT.sumWinit = sum(OUT.subW,2);
mass_removed = (OUT.melt + OUT.moist_sublimation)*1d3;
mass_layer = OUT.subD.*OUT.subZ;
n = 0;
while any(mass_removed>0)
    n = n+1;
    cond1 = mass_removed(:,1)>mass_layer(:,n);
    cond2 = ~cond1 & mass_removed(:,1)>0;
    mass_removed(cond1,1) = mass_removed(cond1,1) - OUT.subD(cond1,n).*OUT.subZ(cond1,n);
    shift_tot(cond1,1) = shift_tot(cond1,1) - OUT.subZ(cond1,n);
    shift_tot(cond2,1) = shift_tot(cond2,1) - mass_removed(cond2,1)./mass_layer(cond2,n).*OUT.subZ(cond2,n);
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
    
    % shift layers up when first layer thickness is completely removed 
    OUT.subZ(i_shift,2:nl-1) = subZ_old(i_shift,3:nl);
    OUT.subT(i_shift,2:nl-1) = subT_old(i_shift,3:nl);
    OUT.subD(i_shift,2:nl-1) = subD_old(i_shift,3:nl);
    OUT.subW(i_shift,2:nl-1) = subW_old(i_shift,3:nl);
    OUT.subZ(i_shift,1) = subZ_old(i_shift,1) + subZ_old(i_shift,2) + shift(i_shift,1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Snow compaction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subD_old = OUT.subD;
subZ_old = OUT.subZ;
mliqmax = zeros(grid.gpsum,grid.nl);


%% Densification of firn [Ligtenberg et al. 2011]
cond = OUT.subD>=C.Dfirn;

grav_const = zeros(grid.gpsum,grid.nl);
grav_const(cond) =  0.07 .* max(1.435-0.151*IN.logyearsnow(cond),0.25) .* (OUT.subD(cond)<550) ...
    + 0.03 .* max(2.366-0.293*IN.logyearsnow(cond),0.25) .* (OUT.subD(cond)>=550);   

OUT.subD(cond) = OUT.subD(cond) + dt./C.yeardays.*grav_const(cond).*IN.yearsnow(cond).*C.g.*(C.Dice-OUT.subD(cond)).* ...
         exp(-C.Ec./(C.rd.*OUT.subT(cond)) + C.Eg./(C.rd.*OUT.subTmean(cond)));
     
OUT.Dens_firn = zeros(grid.gpsum,grid.nl);     
OUT.Dens_firn(cond) = dt./C.yeardays.*grav_const(cond).*IN.yearsnow(cond).*C.g.*(C.Dice-OUT.subD(cond)).* ...
         exp(-C.Ec./(C.rd.*OUT.subT(cond)) + C.Eg./(C.rd.*OUT.subTmean(cond)));
     
     
%% Densification of seasonal snow [Kampenhout et al. 2017]
cond = OUT.subD<C.Dfirn;

% Destructive metamorphism
CC_tap = 200;
CC1 = zeros(grid.gpsum,grid.nl);
CC2 = zeros(grid.gpsum,grid.nl);
CC1(OUT.subD<CC_tap) = 1;
CC1(OUT.subD>=CC_tap) = exp(-0.046.*(OUT.subD(OUT.subD>=CC_tap)-CC_tap));
CC2(OUT.subW==0) = 1;
CC2(OUT.subW~=0) = 2;
CC3 = 2.777d-6;
CC4 = 0.04;

OUT.subD(cond) = OUT.subD(cond) + dt*3600*24.*OUT.subD(cond).*CC3.*CC2(cond).*CC1(cond).*exp(-CC4.*(C.T0-OUT.subT(cond)));
OUT.subD(cond) = min(OUT.subD(cond),C.Dice);

OUT.Dens_destr_metam = zeros(grid.gpsum,grid.nl);   
OUT.Dens_destr_metam(cond) = dt*3600*24.*OUT.subD(cond).*CC3.*CC2(cond).*CC1(cond).*exp(-CC4.*(C.T0-OUT.subT(cond)));

% Overburden pressure
CC5 = 0.1;
CC6 = 0.023;
CC7 = 4.0.*7.62237d6.*OUT.subD./358.0.*1./(1+60.*OUT.subW./(C.Dwater.*OUT.subZ));
Psload = zeros(grid.gpsum,grid.nl);
Visc = zeros(grid.gpsum,grid.nl);
temp = cumsum(OUT.subD.*OUT.subZ,2) - 0.5*OUT.subD.*OUT.subZ;
Psload(cond) = temp(cond); 
Visc(cond) = CC7(cond).*exp(CC5.*(C.T0-OUT.subT(cond))+CC6.*OUT.subD(cond));

OUT.subD(cond) = OUT.subD(cond) + dt*3600*24.*OUT.subD(cond).*Psload(cond)./Visc(cond);
OUT.subD(cond) = min(OUT.subD(cond),C.Dice);

OUT.Dens_overb_pres = zeros(grid.gpsum,grid.nl); 
OUT.Dens_overb_pres(cond) = dt*3600*24.*OUT.subD(cond).*Psload(cond)./Visc(cond);

% Drifting snow
MO = -0.069+0.66.*(1.25-0.0042.*(max(OUT.subD,50)-50));
SI = -2.868.*exp(-0.085.*IN.WS)+1+MO;
cond_SD = SI>0;
z_i = zeros(grid.gpsum,grid.nl);
z_i(:,2:grid.nl) = cumsum(OUT.subZ(:,1:grid.nl-1).*(3.25-SI(:,1:grid.nl-1)),2);
gamma_drift = max(0,SI.*exp(-z_i./0.1));
tau = 48*2*3600;
tau_i = tau./gamma_drift;

OUT.subD(cond & cond_SD) = OUT.subD(cond & cond_SD) + dt*3600*24.*max(350-OUT.subD(cond & cond_SD),0)./tau_i(cond & cond_SD);
OUT.subD(cond & cond_SD) = min(OUT.subD(cond & cond_SD),C.Dice);

OUT.Dens_drift = zeros(grid.gpsum,grid.nl); 
OUT.Dens_drift(cond & cond_SD) = dt*3600*24.*max(350-OUT.subD(cond & cond_SD),0)./tau_i(cond & cond_SD);

%% Update layer thickness, surface height and stored irreducible water after compaction
cond = OUT.subD<C.Dice;

OUT.subZ(cond) = subZ_old(cond).*subD_old(cond)./OUT.subD(cond);

mliqmax(cond) =  OUT.subD(cond).*OUT.subZ(cond).*0.0143.*exp(3.3.*(C.Dice-OUT.subD(cond))./C.Dice)./(1-0.0143.*exp(3.3.*(C.Dice-OUT.subD(cond))./C.Dice))...
    .*0.05.*min(C.Dice-OUT.subD(cond),20);
OUT.subW = min(mliqmax,OUT.subW);

shift(:,1) = sum(OUT.subZ,2)-sum(subZ_old,2);
OUT.surfH(:,1) = OUT.surfH(:,1) + shift(:,1);

runoff_irr = OUT.sumWinit - sum(OUT.subW,2);

%% Update subsurface temperatures after heat conduction / diffusion
dz1 = ((OUT.subZ(:,1)+0.5.*OUT.subZ(:,2)).^2);
dz2 = (0.5.*(OUT.subZ(:,3:grid.nl)+OUT.subZ(:,2:grid.nl-1)).^2);

kk = zeros(grid.gpsum,grid.nl);
c_eff = zeros(grid.gpsum,grid.nl);
kdTdz = zeros(grid.gpsum,grid.nl);

kk = 0.138-1.01d-3.*OUT.subD+3.233d-6.*OUT.subD.^2;
c_eff = OUT.subD .* (152.2+7.122.*OUT.subT); 
    
z_temp = OUT.subZ(:,2:end);
c_eff_temp = c_eff(:,1:end);
kk_temp = kk(:,1:end);
dt_stab = 0.5*min(c_eff_temp,[],2).*min(z_temp,[],2).^2./max(kk_temp,[],2)/3600.0/24.0;

tt = zeros(grid.gpsum,1);
while any(tt<dt)
    subT_old = OUT.subT;
    
    dt_temp = min(dt_stab,-tt+dt);

    tt = tt + dt_temp;
    
    cond_dt = dt_temp>0;
    
    kdTdz(cond_dt,2) = (kk(cond_dt,1).*OUT.subZ(cond_dt,1)+0.5.*kk(cond_dt,2).*OUT.subZ(cond_dt,2)).*(subT_old(cond_dt,2)-OUT.Tsurf(cond_dt,1)) ./ dz1(cond_dt,1);
    kdTdz(cond_dt,3:nl) = (kk(cond_dt,2:nl-1).*OUT.subZ(cond_dt,2:nl-1)+kk(cond_dt,3:nl).*OUT.subZ(cond_dt,3:nl)).*(subT_old(cond_dt,3:nl)-subT_old(cond_dt,2:nl-1)) ./ dz2(cond_dt,:);

    OUT.subT(cond_dt,2) = subT_old(cond_dt,2) + 24.*3600.*dt_temp(cond_dt,1).*(kdTdz(cond_dt,3)-kdTdz(cond_dt,2))./ ...
                        (c_eff(cond_dt,2) .* ...
                        (0.5*OUT.subZ(cond_dt,1)+0.5*OUT.subZ(cond_dt,2)+0.25*OUT.subZ(cond_dt,3)));
    OUT.subT(cond_dt,3:nl-1) = subT_old(cond_dt,3:nl-1) + 24.*3600.*dt_temp(cond_dt,1).*(kdTdz(cond_dt,4:nl)-kdTdz(cond_dt,3:nl-1))./ ...
                        (c_eff(cond_dt,3:nl-1) .* ...
                        (0.25*OUT.subZ(cond_dt,2:nl-2)+0.5*OUT.subZ(cond_dt,3:nl-1)+0.25*OUT.subZ(cond_dt,4:nl)));
                    
    OUT.subT(cond_dt,nl) = subT_old(cond_dt,nl) + 24.*3600.*dt_temp(cond_dt,1).*(C.geothermal_flux-kdTdz(cond_dt,nl))./ ...
                        (c_eff(cond_dt,nl) .* ...
                        (0.25*OUT.subZ(cond_dt,nl-1)+0.75*OUT.subZ(cond_dt,nl)));
end
OUT.subT(:,1) = OUT.Tsurf(:,1) + (OUT.subT(:,2)-OUT.Tsurf(:,1))./(OUT.subZ(:,1)+0.5*OUT.subZ(:,2)).*0.5.*OUT.subZ(:,1);
OUT.subT(OUT.subT>C.T0) = C.T0;
OUT.subCeff = c_eff;
OUT.subK = kk;

%% Refreezing of percolating water and irreducible water storage
avail_W =   OUT.melt*1d3 + ...                                    % melt water
            IN.rain.*1d3 + ...                                % rain fall         
            (OUT.moist_condensation-OUT.moist_evaporation).*1d3;    % condensation or evaporation                                   
avail_W = max(avail_W,0);

subW_old = OUT.subW;
OUT.cpi = 152.2+7.122.*OUT.subT;
c1 = OUT.cpi.*OUT.subD.*OUT.subZ.*(C.T0-OUT.subT)./C.Lm;        % cold content limit on refreezing
c2 = OUT.subZ.*(1-OUT.subD./C.Dice).*C.Dice;                % density limit on refreezing
cond1 = c1 >= c2;
cond2 = ~cond1;
Wlim = zeros(grid.gpsum,nl);                            % maximum potential for refreezing
Wlim(cond1) = c2(cond1);
Wlim(cond2) = c1(cond2);
Wlim = max(Wlim,0);
mliqmax = zeros(grid.gpsum,grid.nl);                    % maximum irreducible water storage
noice = OUT.subD<C.Dice-1;
mliqmax(noice) =  OUT.subD(noice).*OUT.subZ(noice).*0.0143.*exp(3.3.*(C.Dice-OUT.subD(noice))./C.Dice)./(1-0.0143.*exp(3.3.*(C.Dice-OUT.subD(noice))./C.Dice))...
    .*0.05.*min(C.Dice-OUT.subD(noice),20);
Wirr = mliqmax-subW_old;

RP = zeros(grid.gpsum,grid.nl);                         % amount of refreezing of percolating water
leftW = zeros(grid.gpsum,1);
avail_W_loc = zeros(grid.gpsum,1);

z0 = C.perc_depth;
percolation = phys.percolation;

zz = cumsum(OUT.subZ,2)-0.5*OUT.subZ;

% description of water percolation
carrot = zeros(grid.gpsum,nl);      

% prob. density. distribution funtions
if  percolation == 1
    carrot(:,1) = 1;
elseif   percolation == 2      % normal law
    carrot = 2*exp( - zz.^2/2/(z0/3)^2 )/(z0/3)/sqrt(2*pi);
elseif   percolation == 3      % linear
    carrot = 2*(z0 - zz) / z0.^2;
    carrot = max(carrot,0);    % eliminate negative values
elseif   percolation == 4      % uniform
    [~, ind] = min(abs(zz - z0));
    carrot(1:ind) = 1/z0; clear ind
end

% probability function: take into account thickness of each layer
carrot = carrot.*OUT.subZ;

% normalize by sum to not loose or gain water
carrot = carrot./sum(carrot,2);

% multiply by available water mass to distribute it along the profile
carrot = carrot.*avail_W;

for n=1:nl
    avail_W_loc = avail_W_loc + carrot(:,n);
    
    % refreeze percolating water and store irreducible water
    cond1 = avail_W_loc>Wlim(:,n);
    
    RP(cond1,n) = Wlim(cond1,n);        % more water than refreezing limit
    leftW(cond1,1) = (avail_W_loc(cond1,1)-Wlim(cond1,n));
    OUT.subW(cond1,n) = subW_old(cond1,n) + min(leftW(cond1,1),Wirr(cond1,n));
    
    RP(~cond1,n) = avail_W_loc(~cond1); % less water than refreezing limit
    OUT.subW(~cond1,n) = subW_old(~cond1,n);
    
    avail_W_loc = avail_W_loc - RP(:,n) - (OUT.subW(:,n)-subW_old(:,n));
    
    % update temperature and density after refreezing
    OUT.subT(:,n) = OUT.subT(:,n) + C.Lm.*RP(:,n)./(OUT.subD(:,n).*OUT.cpi(:,n).*OUT.subZ(:,n));
    OUT.subD(:,n) = OUT.subD(:,n) + RP(:,n)./OUT.subZ(:,n);
end
avail_W = avail_W_loc;

%% Runoff of slush water
%  - according to linear-reservoir model
%  - recharged by nonrefrozen melt/rain water
%  - discharge depends on slush water amount and runoff time-scale
slush_old = sum(OUT.subS,2);
avail_S = 1.0/(1.0+dt/C.Trunoff).*slush_old;
avail_S(avail_S<1d-25) = 0.0;
runoff_slush = slush_old - avail_S;

%% Storage of slush water
avail_W = avail_W + avail_S;
OUT.subS = zeros(grid.gpsum,grid.nl);
slushspace = zeros(grid.gpsum,grid.nl);

for n=nl:-1:1
    if (any(OUT.subD(:,n)<C.Dice) && any(avail_W>0))
        
        % compute available space and fill with slush water
        % some space is already occupied by irreducible water
        slushspace(:,n) = max(OUT.subZ(:,n).*(1-OUT.subD(:,n)./C.Dice).*C.Dwater - OUT.subW(:,n),0);
        cond1 = avail_W>slushspace(:,n);
        OUT.subS(cond1,n) = slushspace(cond1,n);
        OUT.subS(~cond1,n) = avail_W(~cond1);
        avail_W = avail_W - OUT.subS(:,n);
    end
end
runoff_surface = avail_W;


%% Refreezing of slush water
OUT.cpi = 152.2+7.122.*OUT.subT;                            % specific heat capacity
c1 = OUT.cpi.*OUT.subD.*OUT.subZ.*(C.T0-OUT.subT)./C.Lm;        % cold content limit on refreezing
c2 = OUT.subZ.*(1-OUT.subD./C.Dice).*C.Dice;                % density limit on refreezing
cond1 = c1 >= c2;
cond2 = ~cond1;
Wlim = zeros(grid.gpsum,nl);                            % maximum potential for refreezing
Wlim(cond1) = c2(cond1);
Wlim(cond2) = c1(cond2);
RS = zeros(grid.gpsum,grid.nl);                         % amount of slush water refreezing

for n=nl:-1:1
    if (any(OUT.subS(:,n)>0) && any(OUT.subT(:,n)<C.T0))
        avail_W = OUT.subS(:,n);
        Wlim_loc = Wlim(:,n);

        cond = avail_W>Wlim_loc;
        RS(cond,n) = Wlim_loc(cond);
        RS(~cond,n) = avail_W(~cond);

        % update slush water content, temperature and density
        OUT.subS(:,n) = OUT.subS(:,n) - RS(:,n);
        OUT.subT(:,n) = OUT.subT(:,n) + C.Lm.*RS(:,n)./(OUT.subD(:,n).*OUT.cpi(:,n).*OUT.subZ(:,n));
        OUT.subD(:,n) = OUT.subD(:,n) + RS(:,n)./OUT.subZ(:,n);
    end
end

%% Refreezing of irreducible water
OUT.cpi = 152.2+7.122.*OUT.subT;                            % specific heat capacity
c1 = OUT.cpi.*OUT.subD.*OUT.subZ.*(C.T0-OUT.subT)./C.Lm;    % cold content limit on refreezing
c2 = OUT.subZ.*(1-OUT.subD./C.Dice).*C.Dice;                % density limit on refreezing
cond1 = c1 >= c2;
cond2 = ~cond1;
Wlim = zeros(grid.gpsum,nl);                            % maximum potential for refreezing                        
Wlim(cond1) = c2(cond1);
Wlim(cond2) = c1(cond2);
RI = zeros(grid.gpsum,grid.nl);                         % amount of irreducible water refreezing

for n=nl:-1:1
    if (any(OUT.subW(:,n)>0) && any(OUT.subT(:,n)<C.T0))
        
        % refreeze irreducible water
        avail_W = OUT.subW(:,n);
        cond1 = avail_W>Wlim(:,n);
        RI(cond1,n) = Wlim(cond1,n);
        RI(~cond1,n) = avail_W(~cond1);
        
        % update irreducible water content, tempeature and density
        OUT.subW(:,n) = OUT.subW(:,n) - RI(:,n);
        OUT.subT(:,n) = OUT.subT(:,n) + C.Lm.*RI(:,n)./(OUT.subD(:,n).*OUT.cpi(:,n).*OUT.subZ(:,n));
        OUT.subD(:,n) = OUT.subD(:,n) + RI(:,n)./OUT.subZ(:,n);
    end
end

OUT.refr = 1d-3.*(sum(RP,2) +  sum(RS,2) +  sum(RI,2));
OUT.refr_P = 1d-3.*sum(RP,2);
OUT.refr_S = 1d-3.*sum(RS,2);
OUT.refr_I = 1d-3.*sum(RI,2);
OUT.slushw = sum(OUT.subS,2);
OUT.irrw = sum(OUT.subW,2);

%% Determine volumetric water content
OUT.subWvol = OUT.subW.*1d-3./OUT.subZ;

%% Merge / split layers to increase layer thickness at depths defined in grid.split
if grid.doubledepth
    for n=1:length(grid.split)
        % Merge layers (accumulation case)
        split = grid.split(n);
        cond = OUT.subZ(:,split)<=(2.0^(n-1))*grid.max_subZ & grid.mask_short==1;
        subT_old = OUT.subT;
        subD_old = OUT.subD;
        subW_old = OUT.subW;
        subZ_old = OUT.subZ;
        subS_old = OUT.subS;
        OUT.subZ(cond,split-1) = subZ_old(cond,split-1) + subZ_old(cond,split);
        OUT.subW(cond,split-1) = subW_old(cond,split-1) + subW_old(cond,split);
        OUT.subS(cond,split-1) = subS_old(cond,split-1) + subS_old(cond,split);
        OUT.subD(cond,split-1) = (subZ_old(cond,split-1).*subD_old(cond,split-1) + subZ_old(cond,split).*subD_old(cond,split))./(subZ_old(cond,split-1) + subZ_old(cond,split));
        OUT.subT(cond,split-1) = (subZ_old(cond,split-1).*subT_old(cond,split-1) + subZ_old(cond,split).*subT_old(cond,split))./(subZ_old(cond,split-1) + subZ_old(cond,split));
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
        cond = OUT.subZ(:,split-2)>(2.0^(n-1))*grid.max_subZ & grid.mask_short==1;
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
        runoff_irr_deep(cond,1) = runoff_irr_deep(cond,1) + subW_old(cond,nl);
        runoff_slush(cond,1) = runoff_slush(cond,1) + subS_old(cond,nl);     
    end
end

%% Runoff
OUT.runoff_irr_deep_mean = OUT.runoff_irr_deep_mean .* (1-dt./C.yeardays) + runoff_irr_deep .* dt./C.yeardays;
OUT.runoff = 1d-3.*(runoff_surface + runoff_slush + runoff_irr + OUT.runoff_irr_deep_mean);
OUT.runoff_surf = 1d-3.*runoff_surface;
OUT.runoff_slush = 1d-3.*runoff_slush;
OUT.runoff_irr = 1d-3.*runoff_irr;
OUT.runoff_irr_deep = 1d-3.*OUT.runoff_irr_deep_mean;

end

