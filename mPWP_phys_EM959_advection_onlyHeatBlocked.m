function mPWP_phys(met_input_file, profile_input_file, pwp_output_file,Agrow,Amelt,Aicecover,MoveBackT,MoveBackS)


%  function PWP(met_input_file, profile_input_file, pwp_output_file)

%  MATLAB version 1.4 (8 Oct 2001)

%  This main program is used to drive the pwp subroutine (pwpgo)
%  with a daily cycle of heating, or, an arbitrary (observed) time
%  series of surface flux data.

%  Written by Jim Price, april 27, 1989. Please direct any
%  comments or questions to Jim Price, WHOI, Woods Hole, MA
%  02543, USA, tel. 508-289-2526.

%  Version 2:  changed the gradient Ri relaxation method.
%  Version 3:  clarified the bulk Ri relaxation.
%  Version 4:  fixed an error that disabled background diffusion.

%  Last editted on 20 Jan, 1999 by JFP.

%  Translated for MATLAB 5.X by Peter Lazarevich and Scott Stoermer,
%  GSO/URI, 11 July 2001.


%  -- Option to display output during run. --

diagnostics = 0;

%  -- Initialize the user-defined variables. --

global dt dz days depth dt_save lat g cpw rb rg rkz beta1 beta2 udrag Lf m_kt n_kt Lv Cd
global cpa alp bet rho_a rho_w sigma_sb Patm dw ep
global Pb Pw
global dstab ad


%  dt,			time-step increment (seconds)
%  dz,			depth increment (meters)
%  days,  	the number of days to run
%  depth,		the depth to run
%  ad,          depth below which "advection" occurs (m)
%  dt_save,	time-step increment for saving to file (multiples of dt)
%  lat,  		latitude (degrees)
%  g,				gravity (9.8 m/s^2)
%  cpw,			specific heat of water (4183.3 J/kgC)
%  cpa,         specific heat of air (1005 J/kgC)
%  rb,			critical bulk richardson number (0.65)
%  rg,			critical gradient richardson number (0.25)
%  rkz,			background diffusion (0)
%  beta1, 	longwave extinction coefficient (0.6 m)
%  beta2, 	shortwave extinction coefficient (20 m)
%  Lf,      Latent heat of fusion (3.35e5J/kg)
%  Lv,      Latent heat of vaporisation (2.501e6J/kg)
%  sigma_sb,   Stefan-Boltzmann constant (5.67e-8)
%  m_kt,    Co-efficient for power provided by wind, in Kraus-Turner
%  n_kt,    Co-efficient for power provided by buoyancy, in K-T
%  Cd,      drag coefficient
%  rho_a,   density of air (1.275kg/m^3)
%  rho_w,   density of water (1026kg/m^3)
%  alp,   alpha - thermal expansion coefficient (5.82e-5/C)
%  bet,   beta - salinity contraction coefficient (8e-4)
%  Patm,  atmospheric pressure at the surface (101325Pa)
%  dw,      depth scale of dissipation(10m)
%  ep,     epsilon, ratio of molecular weight of water and dry air (0.622)


dt			= 60;
dz			= 1;
% depth		= 450;
%%%%%%%%%%%%%%% Yixi changed the depth, lat, depth of advection and day (some point in 2019)
    %%%%%%%%%%%%%%% Yixi can't remember why she changed it (19/Feb/2022)
dep_max=load(profile_input_file,'profile');
lat 		= nanmedian(dep_max.profile.lat);
t_ad=gt_sg_sub_filter(dep_max.profile.t,10,0);
% ad=min(dep_max.profile.z(abs(t_ad-nanmax(t_ad))<abs(0.1*nanmax(t_ad))));% Yixi thinks that it's for finding the turning point (10/Oct/2022)
% ad=28;% change it to the initial mix-layer depth, and allow it to change in each time step, 11/Oct/2022
ad=1;% change it to the initial mix-layer depth, and allow it to change in each time step, 11/Oct/2022
dep_max=dep_max.profile.z(end);
depth=(floor(dep_max./dz))*dz;
day_max=load(met_input_file,'met');
days 		= floor(day_max.met.day(end)-day_max.met.day(1));
clearvars dep_max day_max t_ad
    %%%%%%%%%%%%%%% Yixi can't remember why she changed it (19/Feb/2022)
%%%%%%%%%%%%%% Yixi changed the depth, lat, depth of advection and day (some point in 2019)

% dt_save     = 1600; 
dt_save     = 180; % dt_save * dt (dt=120 s) = 3 hr, Yixi changed it to get finer-res figure
g			= 9.8;
cpw			= 4183.3;
cpa         = 1005;
rb			= 0.65;
rg			= 0.25;

rkz			= 0.0001;
beta1 	= 0.6;
beta2 	= 20;
Lf      = 3.35e5;
Lv      = 2.501e6;
sigma_sb= 5.67e-8;
m_kt    = 0.4;
n_kt    = 0.18;
Cd      = 0.001;
rho_a   = 1.275;
Patm    = 101325;
rho_w   = 1026;
alp     = 5.82e-5;
bet     = 8e-4;
dw      = 10;
ep      = 0.62197;


%  -- Initialize additional variables. ---

global f ucon

%  f,			Coriolis parameter
%  ucon,	coefficient of inertial-internal wave dissipation

f = 2*7.29E-5*sin(lat*pi/180);

%  -- Load the air/sea flux data. --

disp(['loading ' met_input_file])

load(met_input_file);

%  -- Interpolate the air/sea flux variables at dt resolution. --

%  Set nmet equal to the number of time increments using a resolution of dt.

nmet 	= days*8.64E4/dt;
time	= met.time(1)+(0:(nmet-1))*dt/8.64E4;

% Check record length

if time(end) > met.time(end)
	time = time(time<met.time(end));
	nmet = length(time);
	disp(['Met input shorter than # of days selected, truncating run to ' num2str(nmet*dt/8.64E4) ' day(s)'])
	pause
end

% Check the time-resolution of the inertial period

if dt > abs(1/10*2*pi/f)
	ans = input('Time step, dt, too large to accurately resolve the inertial period. Is this okay? (y/n)','s');
	if ans == 'n'
		error(['Please restart PWP.m with a new dt <= ' num2str(1/(10*f))])
		%break 
	end
end

qi			= interp1(met.time,met.sw,time);
lw          = (interp1(met.time,met.lw,time));
tair		= interp1(met.time,met.tair,time);
shum        = interp1(met.time,met.shum,time);
UU           = interp1(met.time,met.U,time);
tx			= interp1(met.time,met.tx,time);
ty			= interp1(met.time,met.ty,time);
precip	= interp1(met.time,met.precip,time);

disp('loading complete')

%  -- Load initial t,s profile data. --

disp(['loading ' profile_input_file])

load(profile_input_file);

%  -- Interpolate the profile variables at dz resolution. --

global nz z t s d h_i ml_depth ml_index we iter mr qlat qsens qloss qnet
d=[];
t=[];
s=[];

%  Set nz equal to the number of depth increments + 1 using a resolution of dz.

nz	= 1+depth/dz;
z		= ((0:nz-1)*dz)';
iter = 0;

%Thickness and area of ice is set to zero

h_i = 0;
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Yixi changed this one, so that A can be tuned
%%%%%%%%%%%%%%%%%%%%% 27/Feb/2022
% A = 0;
%%%%%%%%%%%%%%%%%%%%% 27/Feb/2022
%%%%%%%%%%%%%%%%%%%%% Yixi changed this one, so that A can be tuned
%%%%%%%%%%%%%%%%%%%%%
ml_depth = dz;
ml_index = 2;

% Check record length

if z(end) > profile.z(end)
	z = z(z<=profile.z(end));
	nz = length(z);
	disp(['Profile input shorter than depth selected, truncating to ' num2str(z(end)) ' meters'])
	pause
end

% Check the depth-resolution of the profile file

profile_increment = (profile.z(end)-profile.z(1))/(length(profile.z)-1);
if dz < profile_increment/5
	ans = input('Depth increment, dz, is much smaller than profile resolution. Is this okay? (y/n)','s');
	if ans == 'n'
		error(['Please restart PWP.m with a new dz >= ' num2str(profile_increment/5)])
		%break
	end
end

t	= interp1(profile.z,profile.t,z);
s	= interp1(profile.z,profile.s,z);
d	= interp1(profile.z,profile.d,z);

disp('loading complete')


%  -- Initialize additional profile variables at dz resolution. --

global u v absrb
u=[];
v=[];

%  u and v, 	east and north current
%  absrb,			absorbtion fraction

u			= zeros(nz,1);
v			= zeros(nz,1);
absrb = absorb(beta1,beta2);


%  -- Define the variables to be saved. --

pwp_output.dt 			= dt;
pwp_output.dz 			= dz;
pwp_output.lat 			= lat;
pwp_output.z				= z;
pwp_output.time			= [];
pwp_output.t				= [];
pwp_output.s				= [];
pwp_output.d				= [];
pwp_output.u				= [];
pwp_output.v				= [];
pwp_output.hi           = [];
pwp_output.ml           = [];
pwp_output.mlI          = [];
pwp_output.we           = [];
pwp_output.mr           = [];
pwp_output.Pb           = [];
pwp_output.Pw           = [];

pwp_output.qloss=[]; % Yixi 12/2018
pwp_output.qsens=[]; % Yixi 12/2018
pwp_output.qlat=[]; % Yixi 12/2018
pwp_output.qnet=[]; % Yixi 23/01/2019
pwp_output.Ai  = []; % Yixi 27/02/2022
pwp_output.AActive  = []; % Yixi 28/06/2022

%  -- Step through the PWP model. --

disp(['STATUS (out of ' int2str(nmet) ' steps):'])

for m = 1:nmet
    if m==35
        m
    end
    
    if isempty(d)==0
        if sum(imag(d))~=0
            m
            break
            
        end
    end
    
    %%%%%%%%%%%%%%%%%% the main line
	Aactive=pwpgo(qi(m),shum(m),tair(m),lw(m),qi(m),UU(m),tx(m),ty(m),precip(m),m,Agrow,Amelt,Aicecover,MoveBackT,MoveBackS);
    %%%%%%%%%%%%%%%%%% the main line
            
%  -- Store variables. --
	
	if mod(m-1,dt_save) == 0
	pwp_output.time(:,end+1)	= time(m);
	pwp_output.t(:,end+1)			= t;
	pwp_output.s(:,end+1)			= s;
    pwp_output.d(:,end+1)			= d;
    pwp_output.u(:,end+1)			= u;
    pwp_output.v(:,end+1)			= v;
    pwp_output.hi(1,end+1)          = h_i;
    pwp_output.ml(1,end+1)          = ml_depth;
    pwp_output.mlI(1,end+1)         = ml_index;
    pwp_output.we(1,end+1)          = we;
    pwp_output.mr(1,end+1)          = mr;
    pwp_output.Pb(1,end+1)          = Pb;
    pwp_output.Pw(1,end+1)          = Pw;
    
    pwp_output.qloss(1,end+1)       =qloss; % Yixi 12/2019
    pwp_output.qsens(1,end+1)       = qsens; % Yixi 12/2019
    pwp_output.qlat(1,end+1)        = qlat; % Yixi 12/2019
    pwp_output.qnet(1,end+1)        = qnet; % Yixi 23/01/2019
	pwp_output.Ai(1,end+1)          = Aicecover; % Yixi 27/01/2022
	pwp_output.AActive(1,end+1)      = Aactive; % Yixi 28/06/2022

    end

%  -- Save variables. --

	if mod(m, 100) == 0
		disp([int2str(m), ' (' sprintf('%2.1f',100*m/nmet), '%)'])
		eval(['save ' pwp_output_file ' pwp_output'])
	end

end

eval(['save ' pwp_output_file ' pwp_output'])
disp(['Model Run Completed - Selected variables saved in file: ' pwp_output_file])

pwp_output

% ---------------------------------------------------------------------------------------------------------------------


function Aactive=pwpgo(qi,shum,tair,lw,sw,UU,tx,ty,precip,m,Agrow,Amelt,Aicecover,MoveBackT,MoveBackS) % Yixi 08/11/2022
%pwpgo(qi,qloss,emp,tx,ty,m)

%  This subroutine is an implementation of the Price, Weller,
%  Pinkel upper ocean model described in JGR 91, C7 8411-8427
%  July 15, 1986 (PWP). This version was coded and documented by
%  Jim Price in April 1989.

%  Edited on 20 September 1993 by JFP to allow for a critical
%  gradient Richardson number other than 1/4, and to implement a
%  different and a priori better means of achieving convergence
%  of gradient ri mixing (see subroutine stir for details).
%  The major difference is that the revised scheme gives a more
%  smoothly varying mixed layer depth over a diurnal cycle.
%  Edited on 14 December 1998 by JFP to clarify the bulk Ri
%  relaxation.

%  This model also implements an energy budget form of an
%  entrainment parameterization that is very similar to that
%  described in Price, Mooers, Van Leer, JPO 8, 4, 582-599 (and
%  references therein). This part of the model should be
%  treated as developmental only, as there are several features
%  that are arbitrary to this model (i.e., the depth of the ml
%  during times of heating can be the grid interval, dz). To
%  use this parameterization set the bulk Richardson number to
%  zero, and set em1, or em2, or em3 to non-zero. The gradient
%  Richardson number can be zero or not.

global dt dz days depth dt_save lat g cpw rb rg beta1 beta2 udrag Lf h_i
global Cd rho_a rho_w dw alp sigma_sb bet Lv ep cpa Patm
global nz z t s d ml_depth mr ad qsens qlat qloss qnet
global u v absrb Pw Pb t_orig s_orig u_prev v_prev
global f ucon m_kt n_kt ml_index we iter dstab rkz ncp

iter = iter+1;

% %'Relax' values back towards original profile
if m == 1;
    t_orig = t;
    s_orig = s;
end

%%%%%%%%%%%%%%%%%%%%%%% this part is masked as the Worldview shows that
%sea ice was there since day 20/Feb (22/Jul)
% %%%%%%%%%%%%%%%%%%%%%% Yixi added for a zero-A before day 4.6 (when the sea ice starts growing)
% if m<6624
%     A=0;
% end
% %%%%%%%%%%%%%%%%%%%%%% Yixi added for a zero-A before day 4.6 (when the sea ice starts growing)
%%%%%%%%%%%%%%%%%%%%%%% this part is masked as the Worldview shows that
%sea ice was there since day 20/Feb (22/Jul)


%%%%%%%%%%%%%%%%%%%%%%%%%
% if m == 1576800;       %change CDW endpoint to LCDW? (15/08/15 - LCB)
%         t_orig(199:end) = 1.16;
%         s_orig = s_orig;
% end
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% % Advection term, using U as xx m/s, and a dT/dx and a dS/dx term of yy s/m
% 
% ad_i = floor(ad/dz);
% 
% UdTdx=0.08*(-0.01);
% 
% if iter > 1
%     dtemp = UdTdx*(t(ad_i:end)-t_orig(ad_i:end));
%     t(ad_i:end) = t(ad_i:end) + dtemp;
%     ds = UdTdx*(s(ad_i:end)-s_orig(ad_i:end));
%     s(ad_i:end) = s(ad_i:end) + ds;
% end

%%% HOWEVER, Yixi wants to use relaxtion instead of advection, 12/Oct/22
ad_i = floor(ad/dz);

if iter > 1
    dtemp = MoveBackT*(t(ad_i:end)-t_orig(ad_i:end));
    t(ad_i:end) = t(ad_i:end) + dtemp;
    ds = MoveBackS*(s(ad_i:end)-s_orig(ad_i:end));
    s(ad_i:end) = s(ad_i:end) + ds;
end
%%%%%%%%%%%%%%%%%%%%%%%

%Apply background diffusion

diffus

%Set surface temperature
ts = t(1);

%%%%%%%%%%%%%%%%%%%%%%%%% Yixi never knows why we calculate this, instead of
%%%%%%%%%%%%%%%%%%%%%%%%% just using ERA5 values
%Calculate net outgoing longwave radiation
lw = [(0.97*(5.67e-8))*((ts+273.15)^4)] - 0.97*lw; 
%%%%%%%%%%%%%%%%%%%%%%%%% Yixi never knows why we calculate this, instead of
%%%%%%%%%%%%%%%%%%%%%%%%% just using ERA5 values

%Calculate heat fluxes
qsens = ((rho_a*cpa*Cd)*UU.*((ts) - tair));
%Using WMO (2008) definition for ew
ew = 6.112*exp((17.62*(ts))/(243.12+ts)); %in hPa
ew = ew*100;
%Using AOMIP definition of specific humidity. Can't find description of where 0.378 comes from
shum_sat = (ep*ew)/(Patm - (0.378*ew));

qlat = ((rho_a*Lv*Cd)*UU.*[shum_sat - shum]);
qloss = lw + qsens + qlat;
qnet = -qloss + sw;

%Calculate evaporation and emp

evap = qlat/(1000*(Lv));
emp = evap - precip;

% Find the mixed layer index,

ml_depth = z(ml_index);
t(1:ml_index) = mean(t(1:ml_index));
s(1:ml_index) = mean(s(1:ml_index));

%  Apply heat and fresh water fluxes to the mixed layer.

%%%%%%%%%%%% Yixi wants to know when the model gives her complex number and
%%%%%%%%%%%% ruins her day (some point in 2019)
%%%%%%%%%%%%%%%%%%% Yixi added another if session here, make sure that if
%%%%%%%%%%%%%%%%%%% it's ice free, the sea ice will not effect the air-sea
% %%%%%%%%%%%%%%%%%%% interaction (1 Mar, 2022)
% if h_i>0
%     t(1:ml_index) = t(1:ml_index)+[1-Aicecover]*(qi*absrb(1:ml_index))*dt./(dz*d(1:ml_index)*cpw);
%     t(1:ml_index) = mean(t(1:ml_index));
%     d(1:ml_index) = sw_dens0(s(1:ml_index),t(1:ml_index));
%     t(1:ml_index) = t(1:ml_index)+[1-Aicecover]*(-qloss)*dt./(ml_depth*d(1)*cpw);
%     s(1:ml_index) = s(1:ml_index)/(1-(1-Aicecover)*emp*dt/ml_depth);
% else
%     t(1:ml_index) = t(1:ml_index)+1*(qi*absrb(1:ml_index))*dt./(dz*d(1:ml_index)*cpw);
%     t(1:ml_index) = mean(t(1:ml_index));
%     d(1:ml_index) = sw_dens0(s(1:ml_index),t(1:ml_index));
%     t(1:ml_index) = t(1:ml_index)+1*(-qloss)*dt./(ml_depth*d(1)*cpw);
%     s(1:ml_index) = s(1:ml_index)/(1-1*emp*dt/ml_depth);
% end
% %%%%%%%%%%%%%%%%%%% Yixi added another if session here, make sure that if
% %%%%%%%%%%%%%%%%%%% it's ice free, the sea ice will not effect the air-sea
% %%%%%%%%%%%%%%%%%%% interaction (1 Mar, 2022)
% %%%%%%%%%%%% Yixi wants to know when the model gives her complex number and
% %%%%%%%%%%%% ruins her day (some point in 2019)


% Yixi decides that it's not a fake ice cover anymore, but a flux test,
% this is a replacement for the commented lines above (22/Sep/2023)
    t(1:ml_index) = t(1:ml_index)+[1-Aicecover]*(qi*absrb(1:ml_index))*dt./(dz*d(1:ml_index)*cpw);
    t(1:ml_index) = mean(t(1:ml_index));
    d(1:ml_index) = sw_dens0(s(1:ml_index),t(1:ml_index));
    t(1:ml_index) = t(1:ml_index)+[1-Aicecover]*(-qloss)*dt./(ml_depth*d(1)*cpw);
    s(1:ml_index) = s(1:ml_index)/(1-(1-0)*emp*dt/ml_depth);
% Yixi decides that it's not a fake ice cover anymore, but a flux test
% this is a replacement for the commented lines above (22/Sep/2023)

%Comment out lines 390-408 if switching off sea ice model, and activate the
%line below (A=0)
% A = 0;
% t_s = t(1);
% t_fp = sw_fp(s(1),1);
% if t(1) < t_fp; %starts forming sea ice, sets temp to freezing point
%     mr = ([(d(1)/1000)*(cpw/Lf)*(t(1)- t_fp)]);
%     t(1:ml_index) = t_fp;
%     s(1:ml_index) = s(1:ml_index)/(1+(A*mr));
%     h_i = h_i - (mr*ml_depth);
%     A = Agrow;
% elseif h_i > 0;
%     mr =[(d(1)/1000)*(cpw/Lf)*(t(1)- t_fp)];
%     t(1:ml_index) = t_fp;
%     s(1:ml_index) = s(1:ml_index)/(1+(A*mr));
%     h_i = h_i - (mr*ml_depth);
%     A = Amelt;
% else
%     mr = 0;
%     A = 0;
%     h_i = 0;
% end

%%%%%%%%%%%% Yixi made a new session for sea ice forming, which allows the sea
%%%%%%%%%%%% ice to form (reject brine and comsume heat), at the same time, the blocking effect to air-sea
%%%%%%%%%%%% transfer (A is parameter we give, and is constant throughout the model) % 25 Feb 2022
t_s = t(1);
t_fp = sw_fp(s(1),1);
% A=0;
if t(1) < t_fp %starts forming sea ice, sets temp to freezing point
    mr = ([(d(1)/1000)*(cpw/Lf)*(t(1)- t_fp)]);
    t(1:ml_index) = t_fp;
    s(1:ml_index) = s(1:ml_index)/(1+(Agrow*mr));
    h_i = h_i - (mr*ml_depth);
    Athermodynamics = Agrow;
elseif h_i > 0
    mr =[(d(1)/1000)*(cpw/Lf)*(t(1)- t_fp)];
    t(1:ml_index) = t_fp;
    s(1:ml_index) = s(1:ml_index)/(1+(Amelt*mr));
    h_i = h_i - (mr*ml_depth);
    Athermodynamics = Amelt;
else
    mr = 0;
%%%%%%%%%%%%%%%%%%%%% 27/Feb/2022
%%%%%%%%%%%%%%%%%%%%% Yixi changed this one, so that A can be tuned
    Athermodynamics = 0;
%%%%%%%%%%%%%%%%%%%%% 27/Feb/2022
%%%%%%%%%%%%%%%%%%%%% Yixi changed this one, so that A can be tuned
    h_i = 0;
end

%%%%%%%%%%%% Yixi made a new session for sea ice forming, which allows the sea
%%%%%%%%%%%% ice to form (reject brine and comsume heat), at the same time, the blocking effect to air-sea
%%%%%%%%%%%% transfer (A is parameter we give, and is constantthroughout the model) % 25 Feb 2022

%  Absorb solar radiation at depth.

t(ml_index+1:nz) = t(ml_index+1:nz)+(1-Aicecover)*qi*absrb(ml_index+1:nz)*dt./(dz*d(ml_index+1:nz)*cpw);

d = sw_dens0(s,t);


%  Compute the density, and relieve static instability, if it occurs.

remove_si;
ml_depth = z(ml_index);

% Apply wind and freshwater flux mixing (Kraus-Turner type)
% Comment out lines 425-450 if switching off KT mixing.

fw_flux = (s(1)*(((1-0)*(emp))-Athermodynamics*(mr*dz/dt))); %20/Dec, Yixi changed A to Aicethermodynamics and Aicecover, as it captures a more reliable ice forming area
U_Cubed = (((rho_a/rho_w)*(Cd))^(3/2))*UU^(3); %U star cubed

Pw = ((2*m_kt)*exp(-ml_depth/dw)*U_Cubed); % Power for mixing supplied by wind 

Bo = (((g*alp)/(rho_w*cpw))*(qloss-(qi*0.45))) - (g*bet*fw_flux); % Buoyancy fluxes

Pb = (ml_depth/2)*((1+n_kt)*Bo - (1-n_kt)*abs(Bo)); % Power for mixing supplied by buoyancy

we = (Pw+Pb)/(ml_depth*(g*alp*(t(1)-t(ml_index+1)) - g*bet*(s(1) - s(ml_index+1))));

if we >=0
    ml_depth_test = ml_depth + we*dt; %use this line as a test of whether vertical entrainment should continue to be calculated.
    while ml_depth_test > (ml_depth + (dz/2))
        ml_index = ml_index +1;
        ml_depth = z(ml_index);
        Pw = ((2*m_kt)*exp(-ml_depth/dw)*U_Cubed);
        Pb = (ml_depth/2)*((1+n_kt)*Bo - (1-n_kt)*abs(Bo));
        we = (Pw+Pb)/(ml_depth*(g*alp*(t(1)-t(ml_index+1)) - g*bet*(s(1) - s(ml_index+1))));
        ml_depth_test = ml_depth +we*dt;
    end
else
    ml_depth = (Pw /(-Bo)); 
end


%Redefine mixed layer index

ml_index = round(ml_depth/dz)+1;
if ml_index <=1
    ml_index = 2;
end
if isnan(ml_index)==1
    ml_index
end
ml_depth = z(ml_index);

%Mix new values down to bottom of mixed layer.
mix5(ml_index)

if ml_index <= 0
    disp('mixed layer index is less than 1')
end

%  Time step the momentum equation.

%  Rotate the current throughout the water column through an
%  angle equal to inertial rotation for half a time step.

ang = -f*dt/2;

rot(ang);

%  Apply the wind stress to the mixed layer as it now exists.

u_prev = u(1);
v_prev = v(1);

du = (tx/(ml_depth*d(1)))*dt;
dv = (ty/(ml_depth*d(1)))*dt;
u(1:ml_index) = u(1:ml_index)+du;
v(1:ml_index) = v(1:ml_index)+dv;

%  Apply drag to the current (this is a horrible parameterization of
%  inertial-internal wave dispersion).

if ucon > 1E-10
	u = u*(1-dt*ucon);
	v = v*(1-dt*ucon);
end

%  Rotate another half time step.

rot(ang);

%  Finished with the momentum equation for this time step.

%  Do the bulk Richardson number instability form of mixing (as in PWP).

if rb > 1E-5
	bulk_mix(ml_index)
end

ml_depth = z(ml_index);

%  Do the gradient Richardson number instability form of mixing.

if rg > 1e-10
	grad_mix(m);
end

ml_depth = z(ml_index);
Aactive=Aicecover;
% ad=ml_depth;


% -------------------------------------------------------------------------


function bulk_mix(ml_index)

global g rb
global nz z d
global u v

rvc = rb;
ml_index_start = ml_index+1;
for j = ml_index_start:nz
	h 	= z(j);
	dd 	= (d(j)-d(1))/d(1);
	dv 	= (u(j)-u(1))^2+(v(j)-v(1))^2;
	if dv == 0
		rv = Inf;
	else
		rv = g*h*dd/dv;
	end
	if rv > rvc
		break
	else
		mix5(j);
        ml_index = ml_index+1;
	end
end


% -------------------------------------------------------------------------


function grad_mix(m)

%  This function performs the gradeint Richardson Number relaxation
%  by mixing adjacent cells just enough to bring them to a new
%  Richardson Number.

    if m==7638
        m
    end


    
global dz g rg
global nz z t s d 
global u v

rc 	= rg;

%  Compute the gradeint Richardson Number, taking care to avoid dividing by
%  zero in the mixed layer.  The numerical values of the minimum allowable
%  density and velocity differences are entirely arbitrary, and should not
%  effect the calculations (except that on some occasions they evidnetly have!)

j1 = 1;
j2 = nz-1;

count_loop=0;
if m==7638
m
end
while 1 % loop forever , Yixi commented, 2nd Jun, 2019 (Yixi was very desperate 
    % working in a super hot Sunday debugging the model - as the model gave 
    % her complex which totally pissed her off)
    count_loop=count_loop+1;
   
    if m==7638
        if count_loop==2
            count_loop
        end
        if count_loop==3
            count_loop
        end
        if count_loop==4
            count_loop
        end
        if count_loop==5
            count_loop
        end
        if count_loop==6
            count_loop
        end
        if count_loop==7
            count_loop
        end
        if count_loop==8
            count_loop
        end
        if count_loop==9
            count_loop
        end
        if count_loop==10
            count_loop
        end
    end
        
	for j = j1:j2
        
%         if j==20
%             j
%         end
        
            if sum(imag(d))~=0
                j                
            end
        
		if j <= 0
			keyboard
		end
		dd = (d(j+1)-d(j))/d(j);
		dv = (u(j+1)-u(j))^2+(v(j+1)-v(j))^2;
		if dv < 1e-10
			r(j) = Inf;
		else
			r(j) = g*dz*dd/dv;
		end
	end

	%  Find the smallest value of r in profile

	rs = min(r);
	js = min(find(r==rs));
    
	%  Check to see whether the smallest r is critical or not.
    
	if rs > rc
		return
	end

	%  Mix the cells js and js+1 that had the smallest Richardson Number

	stir(rc,rs,js);
    remove_si 

	%  Recompute the Richardson Number over the part of the profile that has changed

	j1 = js-2;
	if j1 < 1
		 j1 = 1;
	end
	j2 = js+2;
	if j2 > nz-1
		 j2 = nz-1;
	end
end


% -------------------------------------------------------------------------


function a = stir(rc,r,j)

%  This subroutine mixes cells j and j+1 just enough so that
%  the Richardson number after the mixing is brought up to
%  the value rnew. In order to have this mixing process
%  converge, rnew must exceed the critical value of the
%  richardson number where mixing is presumed to start. If
%  r critical = rc = 0.25 (the nominal value), and r = 0.20, then
%  rnew = 0.3 would be reasonable. If r were smaller, then a
%  larger value of rnew - rc is used to hasten convergence.

%  This subroutine was modified by JFP in Sep 93 to allow for an
%  aribtrary rc and to achieve faster convergence.

global t s d u v

rcon 			= 0.02+(rc-r)/2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rnew 			= rc+rcon/5;
f 				= 1-r./rnew;
dt				= (t(j+1)-t(j))*f/2;
t(j+1)		= t(j+1)-dt;
t(j) 			= t(j)+dt;
ds				= (s(j+1)-s(j))*f/2;
s(j+1)		= s(j+1)-ds;
s(j) 			= s(j)+ds;
d(j:j+1)	= sw_dens0(s(j:j+1),t(j:j+1));
du				= (u(j+1)-u(j))*f/2;
u(j+1)		= u(j+1)-du;
u(j) 			= u(j)+du;
dv				= (v(j+1)-v(j))*f/2;
v(j+1)		= v(j+1)-dv;
v(j) 			= v(j)+dv;



% -------------------------------------------------------------------------


function mix5(j)

%  This subroutine mixes the arrays t, s, u, v down to level j.

global t s d u v 

t(1:j) = mean(t(1:j));
s(1:j) = mean(s(1:j));
d(1:j) = sw_dens0(s(1:j),t(1:j));
u(1:j) = mean(u(1:j));
v(1:j) = mean(v(1:j));



% -------------------------------------------------------------------------


function rot(ang)

%  This subroutine rotates the vector (u,v) through an angle, ang

global u v

r = (u+i*v)*exp(i*ang);
u = real(r);
v = imag(r);


% -------------------------------------------------------------------------


function remove_si

%  Find and relieve static instability that may occur in the
%  density array d. This simulates free convection.
%  ml_index is the index of the depth of the surface mixed layer after adjustment,

global d ml_index

% while 1
% 	ml_index = min(find(diff(d)<0));
% 	if isempty(ml_index)
% 		break
% 	end
% 	mix5(ml_index+1);
% end

while d(ml_index) > d(ml_index+1)
    mix5(ml_index+1)
    ml_index = ml_index+1;
end

% -------------------------------------------------------------------------


function absrb = absorb(beta1,beta2)

%  Compute solar radiation absorption profile. This
%  subroutine assumes two wavelengths, and a double
%  exponential depth dependence for absorption.

%  Subscript 1 is for red, non-penetrating light, and
%  2 is for blue, penetrating light. rs1 is the fraction
%  assumed to be red.

global nz dz

rs1 = 0.6;
rs2 = 1.0-rs1;
absrb = zeros(nz,1);
z1 = (0:nz-1)*dz;
z2 = z1 + dz;
z1b1 = z1/beta1;
z2b1 = z2/beta1;
z1b2 = z1/beta2;
z2b2 = z2/beta2;
absrb = (rs1*(exp(-z1b1)-exp(-z2b1))+rs2*(exp(-z1b2)-exp(-z2b2)))';
abcde = 1;

% -------------------------------------------------------------------------

% 
% function a = diffus(dstab,a)
% 
% %  This subroutine applies a simple diffusion
% %  operation to the array a. It leaves the endpoints
% %  unchanged (assumes nothing about the
% %  boundary conditions).
% 
% global nz
% 
% a(2:nz-1) = a(2:nz-1)+dstab*(a(1:nz-2)-2*a(2:nz-1)+a(3:nz));

%-------------------------------------------------------------------------

function diffus

global t s u v d nz
%dstab = (dt/dz^2)*K
dstab1 = 0.001;
dstab2 = 0.005;
dstab3 = 0.001;
x = 1;

while x == 1;
t(2:nz-1) = t(2:nz-1)+dstab1*(t(1:nz-2)-2*t(2:nz-1)+t(3:nz));
s(2:nz-1) = s(2:nz-1)+dstab1*(s(1:nz-2)-2*s(2:nz-1)+s(3:nz));
d = sw_dens0(s,t);
u(2:nz-1) = u(2:nz-1)+dstab2*(u(1:nz-2)-2*u(2:nz-1)+u(3:nz));
v(2:nz-1) = v(2:nz-1)+dstab2*(v(1:nz-2)-2*v(2:nz-1)+v(3:nz));
x = x+1;
end

% -------------------------------------------------------------------------


function dens = sw_dens0(S,T)

% SW_DENS0   Denisty of sea water at atmospheric pressure
%=========================================================================
% SW_DENS0  $Revision: 1.3 $  $Date: 1994/10/10 04:54:09 $
%           Copyright (C) CSIRO, Phil Morgan 1992
%
% USAGE:  dens0 = sw_dens0(S,T)
%
% DESCRIPTION:
%    Density of Sea Water at atmospheric pressure using
%    UNESCO 1983 (EOS 1980) polynomial.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (IPTS-68)]
%
% OUTPUT:
%   dens0 = density  [kg/m^3] of salt water with properties S,T,
%           P=0 (0 db gauge pressure)
% 
% AUTHOR:  Phil Morgan 92-11-05  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%     Unesco 1983. Algorithms for computation of fundamental properties of 
%     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%
%     Millero, F.J. and  Poisson, A.
%     International one-atmosphere equation of state of seawater.
%     Deep-Sea Res. 1981. Vol28A(6) pp625-629.
%=========================================================================

% CALLER: general purpose, sw_dens.m
% CALLEE: sw_smow.m

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=2
   error('sw_dens0.m: Must pass 2 parameters')
end %if

[mS,nS] = size(S);
[mT,nT] = size(T);

if (mS~=mT) | (nS~=nT)
   error('sw_dens0.m: S,T inputs must have the same dimensions')
end %if

Transpose = 0;
if mS == 1  % a row vector
  S = S(:);
  T = T(:);
  Transpose = 1;
end %if

%----------------------
% DEFINE CONSTANTS
%----------------------
%     UNESCO 1983 eqn(13) p17.

b0 =  8.24493e-1;
b1 = -4.0899e-3;
b2 =  7.6438e-5;
b3 = -8.2467e-7;
b4 =  5.3875e-9;

c0 = -5.72466e-3;
c1 = +1.0227e-4;
c2 = -1.6546e-6;

d0 = 4.8314e-4;

%$$$ dens = sw_smow(T) + (b0 + b1*T + b2*T.^2 + b3*T.^3 + b4*T.^4).*S  ...
%$$$                    + (c0 + c1*T + c2*T.^2).*S.*sqrt(S) + d0*S.^2;
if S<0
    S
end
dens = sw_smow(T) + (b0 + (b1 + (b2 + (b3 + b4*T).*T).*T).*T).*S  ...
                   + (c0 + (c1 + c2*T).*T).*S.*sqrt(S) + d0*S.^2;	% complex - Yixi wants to know where is the complex number that fuck her life       

if Transpose
  dens = dens';
end %if

return


% -------------------------------------------------------------------------


function dens = sw_smow(T)

% SW_SMOW    Denisty of standard mean ocean water (pure water)
%=========================================================================
% SW_SMOW  $Revision: 1.3 $  $Date: 1994/10/10 05:51:46 $
%          Copyright (C) CSIRO, Phil Morgan 1992.
%
% USAGE:  dens = sw_smow(T)
%
% DESCRIPTION:
%    Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980. 
%
% INPUT: 
%   T = temperature [degree C (IPTS-68)]
%
% OUTPUT:
%   dens = density  [kg/m^3] 
% 
% AUTHOR:  Phil Morgan 92-11-05  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%     Unesco 1983. Algorithms for computation of fundamental properties of 
%     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%     UNESCO 1983 p17  Eqn(14)
%
%     Millero, F.J & Poisson, A.
%     INternational one-atmosphere equation of state for seawater.
%     Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
%=========================================================================

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
% TEST INPUTS
if nargin ~= 1
   error('sw_smow.m: Only one input argument allowed')
end %if

Transpose = 0;
[mT,nT] = size(T);
if mT == 1 % a row vector
   T = T(:);
   Tranpose = 1;
end %if

%----------------------
% DEFINE CONSTANTS
%----------------------
a0 = 999.842594;
a1 =   6.793952e-2;
a2 =  -9.095290e-3;
a3 =   1.001685e-4;
a4 =  -1.120083e-6;
a5 =   6.536332e-9;

dens = a0 + (a1 + (a2 + (a3 + (a4 + a5*T).*T).*T).*T).*T;

if Transpose
  dens = dens';
end %if

return

%------------------------------------------------------------------------


function c = sw_satO2(S,T)

% SW_SATO2   Satuaration of O2 in sea water
%=========================================================================
% sw_satO2 $Id: sw_satO2.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
%          Copyright (C) CSIRO, Phil Morgan 1998.
%
% USAGE:  satO2 = sw_satO2(S,T)
%
% DESCRIPTION:
%    Solubility (satuaration) of Oxygen (O2) in sea water
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (ITS-68)]
%
% OUTPUT:
%   satO2 = solubility of O2  [ml/l]
%
% AUTHOR:  Phil Morgan 97-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)
%
%$$$ #include "disclaimer_in_code.inc"
%
% REFERENCES:
%    Weiss, R. F. 1970
%    "The solubility of nitrogen, oxygen and argon in water and seawater."
%    Deap-Sea Research., 1970, Vol 17, pp721-735.
%=========================================================================

% Modifications
% 99-06-25. Lindsay Pender, Fixed transpose of row vectors.
% 03-12-12. Lindsay Pender, Converted to ITS-90.

% CALLER: general purpose
% CALLEE:

%$$$ #ifdef VARIANT_PRIVATE
%$$$ %***********************************************************
%$$$ %$Id: sw_satO2.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
%$$$ %
%$$$ %$Log: sw_satO2.m,v $
%$$$ %Revision 1.1  2003/12/12 04:23:22  pen078
%$$$ %*** empty log message ***
%$$$ %

%$$$ %
%$$$ %***********************************************************
%$$$ #endif

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=2
   error('sw_satO2.m: Must pass 2 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);


% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
   error('sw_satO2: S & T must have same dimensions')
end %if

%------
% BEGIN
%------

% convert T to Kelvin
T = 273.15 + T * 1.00024;

% constants for Eqn (4) of Weiss 1970
a1 = -173.4292;
a2 =  249.6339;
a3 =  143.3483;
a4 =  -21.8492;
b1 =   -0.033096;
b2 =    0.014259;
b3 =   -0.0017000;

% Eqn (4) of Weiss 1970
lnC = a1 + a2.*(100./T) + a3.*log(T./100) + a4.*(T./100) + ...
      S.*( b1 + b2.*(T./100) + b3.*((T./100).^2) );

c = exp(lnC);


