%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% UNATTENUATED INCOMING SOLAR RADIATION & SHADING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT] = TIME_flux_insol(grid,time,OUT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Top-of-atmosphere radiation, solar declination and hour angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time as a fraction of a year (in radians)
trad = 2*pi*((time.TCUR-datenum(year(datetime(datevec(time.TCUR))),1,1))...
    /365.242);
    
% Top-of-atmosphere radiation on a surface normal to incident beam
OUT.I0 = 1353.0*(1.0+0.034*cos(trad));                                      % SOURCE: Meyers and Dale (1983)                                          
    
% Solar declination (d)
d = 0.322003-22.971*cos(trad)-0.357898*cos(2d0*trad) ...                    % SOURCE: Iqbal (1983)
            - 0.14398*cos(3d0*trad)+3.94638*sin(trad) ...
			+ 0.019334*sin(2d0*trad)+0.05928*sin(3d0*trad);


% Solar hour angle (h)
B = 360/365.*(time.TCUR-datenum(year(datetime(datevec(time.TCUR))),1,1)...
    -81);
Tcor_ecc = 9.87.*sin(2.*B*pi/180)-7.53.*cos(B*pi/180)-1.5*sin(B*pi/180);    % correction for eccentricity 
Tcor_lon = 4*(grid.lon-15*time.dT_UTC);                                     % correction for longitude within time-zone
Tcor = Tcor_ecc + Tcor_lon;
LST =   hour(datetime(datevec(time.TCUR))) + ...
        minute(datetime(datevec(time.TCUR))) / 60 + ...
        Tcor / 60;                       
h = 15*(LST-12);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Incoming solar radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On a flat surface
OUT.TOAflat = OUT.I0.*(sin(grid.lat*pi/180.0).*sin(d*pi/180.0) + ...        % SOURCE: Iqbal (1983)
    cos(grid.lat*pi/180.0).*cos(d*pi/180.0).*cos(h*pi/180.0));

% On an inclined surface
slopebeta = grid.slope_beta;
slopegamma = grid.slope_gamma; 
OUT.TOA = OUT.I0.*((sin(grid.lat*pi/180.0).*cos(slopebeta) ...        
        - cos(grid.lat*pi/180.0).*sin(slopebeta).*cos(slopegamma)) ...
        .*sin(d*pi/180.0)+ (cos(grid.lat*pi/180.0).*cos(slopebeta) ...
        + sin(grid.lat*pi/180.0).*sin(slopebeta).*cos(slopegamma)) ...
        .*cos(d*pi/180.0).*cos(h*pi/180.0) + cos(d*pi/180.0).*sin( ...
        slopebeta).*sin(slopegamma).*sin(h*pi/180.0));
OUT.TOA(OUT.TOA<0) = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Shading by the surrounding topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elevationangle = asin(sin(grid.lat.*pi/180.0).*sin(d*pi/180d0)...           % SOURCE: Iqbal (1983)
    +cos(grid.lat.*pi/180.0).*cos(d*pi/180.0).*cos(h*pi/180.0));

if h<0
    azimuth = acos((cos(h*pi/180.0).*cos(d*pi/180.0).*sin(grid.lat*... 
        pi/180.0) - sin(d*pi/180.0).*cos(grid.lat*pi/180.0))./cos(...
        elevationangle));
else
    azimuth = -acos((cos(h*pi/180.0).*cos(d*pi/180.0).*sin(grid.lat... 
        *pi/180.0)- sin(d*pi/180.0).*cos(grid.lat*pi/180.0))./cos(...
        elevationangle));
end
az = real(azimuth);

yl = grid.Ly;
xl = grid.Lx;

ddx(az<=-0.75*pi) = -tan(pi+az(az<=-0.75*pi));
ddx(az<=-0.25*pi & az>-0.75*pi) = -1;
ddx(az<=0.25*pi & az>-0.25*pi) = tan(az(az<=0.25*pi & az>-0.25*pi));
ddx(az<=0.75*pi & az>0.25*pi) = 1;
ddx(az>0.75*pi) = tan(pi-az(az>0.75*pi));

ddy(az<=-0.75*pi) = 1;
ddy(az<=-0.25*pi & az>-0.75*pi) = -tan(0.5*pi+az(az<=-0.25*pi & az>...
    -0.75*pi));
ddy(az<=0.25*pi & az>-0.25*pi) = -1;
ddy(az<=0.75*pi & az>0.25*pi) = -tan(0.5*pi-az(az<=0.75*pi & az>...
    0.25*pi));
ddy(az>0.75*pi) = 1;

azmean = nanmean(az(:));

istart = 1*(azmean<=0) + yl*(azmean>0);
jstart = xl*(azmean<=-0.5*pi || azmean>0.5*pi) + 1*(azmean>-0.5*pi ...
    && azmean<=0.5*pi);
iend = yl*(azmean<=0) + 1*(azmean>0);
jend = 1*(azmean<=-0.5*pi || azmean>0.5*pi) + xl*(azmean>-0.5*pi ...
    && azmean<=0.5*pi);
ii = 1*(azmean<=0) - 1*(azmean>0);
jj = -1*(azmean<=-0.5*pi || azmean>0.5*pi) + 1*(azmean>-0.5*pi ...
    && azmean<=0.5*pi);

ddx_2D = zeros(xl,yl);
ddy_2D = zeros(xl,yl);
elevationangle_2D = zeros(xl,yl);
for n=1:length(grid.xind)
    ddx_2D(grid.xind(n),grid.yind(n)) = ddx(n);
    ddy_2D(grid.xind(n),grid.yind(n)) = ddy(n);
    elevationangle_2D(grid.xind(n),grid.yind(n)) = elevationangle(n);
end

shade_2D = ones(xl,yl);
for i=istart:ii:iend
    for j=jstart:jj:jend
        if (grid.mask_2D(j,i)==1)
            count=1;
            kk=round(i+ddx_2D(j,i)*count);
            ll=round(j+ddy_2D(j,i)*count);
            while (kk~=0 && kk~=yl+1 && ll~=0 && ll~=xl+1)
                griddist = sqrt((grid.x_2D(ll,kk)-grid.x_2D(j,i))^2 + ...
                    (grid.y_2D(ll,kk)-grid.y_2D(j,i))^2);
                gridangle_2D = atan((grid.z_2D(ll,kk)-grid.z_2D(j,i))/griddist);
                if (elevationangle_2D(j,i)<=gridangle_2D)
                    shade_2D(j,i)=1;
                    kk=0;
                elseif (shade_2D(ll,kk)==0)
                    shade_2D(j,i)=0;
                    kk=0;
                else
                    shade_2D(j,i)=0;
                    count=count+1;
                    kk=round(i+ddx_2D(j,i)*count);
                    ll=round(j+ddy_2D(j,i)*count);
                end
                if(elevationangle_2D(j,i)<0),shade_2D(j,i)=1; end
            end
        end
    end
end
shade_2D(:,1) = 0;
shade_2D(1,:) = 0;
shade_2D(:,grid.Ly) = 0;
shade_2D(grid.Lx,:) = 0;
OUT.shade = reshape(shade_2D,[grid.Lx*grid.Ly,1]);
OUT.shade = OUT.shade(grid.mask_2D==1);


end
