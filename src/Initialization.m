 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global COL ROW
global INDC1 INDC2
global tpw tvisw npvtw
global tpg tvisg npvtg  
global dpth vb dx dy dz
global k_temp kx 
global krw krg
global cnstx cnstTx 
global pw pw_0 pcgw pcgw_0 Pd0
global phi_0 phi
global sw sw_0 sh sh_0 sg sg_0 
global cl cl_0 cm cm_0
global dnw dnw_0 dng dng_0 dnh dns dni
global vsw vsw_0 vsg vsg_0
global T T_0
global Mh2o Mch4 Msalt
global lambda lambda_g lambda_h lambda_w lambda_s 
global Cp_g Cp_h Cp_w Cp_s 
global  swr wn gn
global eff_stress 
global qw qg qt
global surf_T  Tgrad seafloor
global g  
global Dmethanex  Dsaltx 
global sgr
global total_stress
global p_salinity T_salinity salinity
global Fracture w0 C_frac S_frac Pcr 
global Heat_Capacity_ini kx_ini phi_ini pcgw_ini 
global Cl_Bottom
global hydrate_on 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.81;     % gravitational constant
dnh=912;    % density of hydrate
dns=2700;   % density of solid grain
dni=917;    % ice density at 0 degree c
lithostat=1600; % bulk density

Mh2o=5.75*0.018/(0.016+5.75*0.018); % molar weight fraction of water in hydrate
Mch4=0.016/(0.016+5.75*0.018);      % molar weight fraction of methane in hydrate
Msalt=0.05844;                      % molecular weight of salt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Thermal Constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cp_g=2170;  % KJ/K*kg - heat capacity of methane gas
Cp_h=2200;  % KJ/K*kg - heat capacity of hydrate
Cp_w=4.2e3; % KJ/K*kg - heat capacity of water
Cp_s=1.381e3; % KJ/K*kg - heat capacity of sediment grains

lambda_g=0.034; % KW/m*K - thermal conductivity of gas
lambda_h=0.60;  % KW/m*K - thermal conductivity of hydrate
lambda_w=0.56;  % KW/m*K - thermal conductivity of water
lambda_s=1.6;   % KW/m*K - thermal conductivity of sediment grains

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Phase Boundary for Methane Hydrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up P,T,S diagram
load p_salinity_8_8_2023.mat
T_salinity=(-20:0.01:40);
salinity=(0:0.5:22)/100;

hydrate_on = 0;     % no hydrate = 0; with hydrate = 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%        Viscosity of Water and Gas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npvtw=10;
tpw=[14.7 270 520 1015 2015 2515 3015 4015 5015 9015]*6895;
tvisw=[1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31 1.31]*1E-3;

npvtg=10;
tpg=[14.7 270 520 1015 2015 2515 3015 4015 5015 9015]*6895;
tvisg=ones(1,10)*0.02*1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Relative Permeability 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sgr=0;%0.02;
swr=0;
wn=3.8;   % constant used in equations to calculate relative water permeability 
gn=2.1;   % constant used in equations to calculate relative gas permeability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Fracture Properties 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w0=50e-6;                       % maximum fracture size in m
S_frac=1000e-6;                 % fracture spacing in m 
C_frac=1/(1.8e10*4.65e-5)*0.5;  % fracture compressibility in 1/Pa     
Pcr=1e6;                        % fracture tensile strength

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seafloor=1000;          % Water depth
model_top=seafloor;     % Top of model 

ROW=1;
COL=56;                 % Total number of grids
dx=20.*ones(COL,ROW);   % Vertical grid size
dy=10*ones(COL,ROW);
dz=10*ones(COL,ROW);

% Vertical depth from sea level 
dpth=zeros(COL,1);
dpth(1)=model_top;
for i=1:COL-1
    dpth(i+1)=dpth(i)+(dx(i,1)+dx(i+1,1))*0.5;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Boundary Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qw=1;               % Bottom water flux in kg/m2/yr 
qg=0;               % Bottom gas flux in kg/m2/yr 
qt=0.0474;          % Geothermal heat flux in W/m2
Cl_Bottom=0.16;     % Salinity of brine injected from bottom boundary in wt.%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Initialial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------initial saturations--------------------------------
sw_0=ones(COL,ROW);     % initial water saturation
sh_0=zeros(COL,ROW);    % initial hydrate saturation
sg_0=zeros(COL,ROW);    % initial gas saturation

sw=sw_0;
sh=sh_0;
sg=sg_0;

%-------------initial porosity, permeability, capillary entry pressure ----

phi_mud = 0.4;      % sediment matrix porosity
perm_mud=1e-19;     % sediment matrix permeability in m-2
pcgw_mud = 1e6;     % sediment matrix capillary entry pressure
mud_bottom = 51;    % the grid number at the base of mud

phi_0 = phi_mud.*ones(COL,ROW);     % initial sediment porosity
k_temp = perm_mud.*ones(COL,ROW);   % initial sediment permeability
Pd0 = pcgw_mud.*ones(COL,ROW);      % initial sediment capillary entry pressure
kx = k_temp;
pcgw = Pd0; 

kx_ini = kx;
phi_ini=phi_0;
pcgw_ini=pcgw;

% ------- Harmonic averaging of permeability and flow operator-------------

for i=1:COL
    for j=1:ROW
        vb(i,j)=dx(i,j)*dy(i,j)*dz(i,j);    % Grid volume
        if i < COL
            cnstx(i,j)=2.0*dz(i+1,j)*dy(i+1,j)*dz(i,j)*dy(i,j)*kx(i+1,j)*kx(i,j)/ ...
                (dx(i,j)*dz(i+1,j)*dy(i+1,j)*kx(i+1,j) + dx(i+1,j)*dz(i,j)*dy(i,j)*kx(i,j)); 
            cnstTx(i,j)=2.0*dz(i+1,j)*dy(i+1,j)*dz(i,j)*dy(i,j)/ ...
                (dx(i,j)*dz(i+1,j)*dy(i+1,j) + dx(i+1,j)*dz(i,j)*dy(i,j));
        else
            cnstx(i,j)=0.0;
            cnstTx(i,j)=0.0;
        end
    end
end

%-----------------------initial temeprature--------------------------------

surf_T=5;   % seafloor temperature in degC
Tgrad=0.04; % Initial temperature gradient in degC/m
T_0=surf_T+(dpth-seafloor).*Tgrad; % Initial  temperature
T = T_0;

%--------------------------initial salinity-------------------------------

cl_0=0.035.*ones(COL,ROW);
cl=cl_0;

%----------------------initial pressure, effective and total streess-------
for i=1:COL
    for j=1:ROW
        pw_0(i,j)=dpth(i,j)*1033.6*g;   %initial hydrastatic pressure
        eff_stress(i,j)=pw_0(1,j)+lithostat*g*(dpth(i,j)-dpth(1,j))-pw_0(i,j); % initial effective stress       
    end
end
pw = pw_0;
eff_stress(1)=0;
total_stress=pw_0+eff_stress;

%------------------------- initial fracture distribution ------------------

Fracture=zeros(COL,ROW);
 
%-------- initial methane concentration, density, viscosity and porosity---
for i=1:COL
    for j=1:ROW
         if sh_0(i,j)>0
            cm_0(i,j)=hydrate_solubility(pw_0(i,j)/1e6,T_0(i,j),cl_0(i,j)/Msalt/(1-cl_0(i,j)))*0.016;
        elseif sg_0(i,j)>0
            cm_0(i,j)=methane_solubility(pw_0(i,j)/1e6,T_0(i,j),cl_0(i,j)/Msalt/(1-cl_0(i,j)))*0.016;
         else
             cm_0(i,j)=0;
        end
        
        dnw_0(i,j)=brine_density(pw_0(i,j)/1e6,T_0(i,j),cl_0(i,j)/Msalt/(1-cl_0(i,j)), cm_0(i,j)/0.016);       % initial brine water density
        vsw_0(i,j)=interp(tpw,tvisw,npvtw,pw_0(i,j),1); % water viscosity
        dng_0(i,j)=gas_density(pw_0(i,j)/1e6,T_0(i,j));   % initial gas density
        vsg_0(i,j)=interp(tpg,tvisg,npvtg,pw_0(i,j),1); % gas viscosity
    end
end
cm=cm_0;
dnw=dnw_0;
vsw=vsw_0;
dng=dng_0;
vsg=vsg_0;

%------------ update total stress & effective stress ------------

total_stress=dnw_0(1,1)*g*seafloor.*ones(COL,ROW);
for i=2:COL
    for j=1:ROW
        total_stress(i,j)=total_stress(i-1,j)+(phi_0(i,j)*dnw_0(i,j)+(1-phi_0(i,j))*dns)*g*(dpth(i,j)-dpth(i-1,j));
    end
end
eff_stress=total_stress-pw_0;
eff_stress(1) = 0;

%------------ heat capacity and conductivity of bulk system ---------------

for i=1:COL
    for j=1:ROW
        lambda(i,j)=(1-phi_0(i,j))*lambda_s+phi_0(i,j)*(sw_0(i,j)*lambda_w+sg_0(i,j)*lambda_g+sh_0(i,j)*lambda_h);
    end
end
Heat_Capacity_ini=phi_0.*dnw_0*Cp_w + (1-phi_0)*dns*Cp_s;

%------------ initial porosity for flow: fracture porosity ----------------

phi_0(1:mud_bottom,:) = 2*w0 /S_frac.*ones(mud_bottom,1);
phi=phi_0;

%-----------------initial water & gas relative permeability----------------

for i=1:COL
    for j=1:ROW
        if sg(i,j)/(1-sh(i,j))>sgr
            krg(i,j)=(sg(i,j)/(1-sh(i,j))-sgr)^gn;     % initial gas phase relative permeability
        else
            krg(i,j)=0;
        end
        if sw(i,j)/(1-sh(i,j))>swr
            krw(i,j)=(sw(i,j)/(1-sh(i,j))-swr)^wn;      % initial water phase relative permeability
        else
            krw(i,j)=0;
        end
    end
end

% ------------------initial index (INDC1, INDC2) values--------------------

for i=1:COL
    for j=1:ROW
        % INDC1=1, hydrate stable zone; INDC1=0, gas stable zone; phase
        % boundary, INDC1=0.
       
        % INDC2=1, one phase present; INDC2=3, two phases present; INDC2=5,
        % three phases present;
        
         if cl(i,j) <interpolation2(p_salinity,T_salinity,salinity,pw(i,j)/1e6,T(i,j))
            INDC1(i,j)=1;
        else
            INDC1(i,j)=0;
        end
        
        INDC2(i,j)=1;
        if sg(i,j)>0 || sh(i,j)>0
            INDC2(i,j)=3;
        elseif sg(i,j)>0 && sh(i,j)>0
            INDC2(i,j)=5;
        end
    end
end
%----------------initial thermal conductivity, gas and water mobility-----------------------------

for i=1:COL
    for j=1:ROW
        calc_trans(i,j);    
    end
end

%---------------------- diffusion coefficient ----------------------------

D0=3.6e-10.*ones(COL,ROW);  % Efective diffusion coefficient 
Dmethanex=zeros(COL,ROW);   % Efective diffusion coefficient of dissolved methane
Dsaltx=zeros(COL,ROW);      % Efective diffusion coefficient of salt
for i=1:COL-1
    Dmethanex(i,j)=dy(i,j)*dz(i,j)/((dx(i,j)+dx(i+1,j))*0.5)*3.6e-10;
    Dsaltx(i,j)=dy(i,j)*dz(i,j)/((dx(i,j)+dx(i+1,j))*0.5)*3.6e-10;
end
Dmethanex(COL,:)=Dmethanex(COL-1,:);
Dsaltx(COL,:)=Dsaltx(COL-1,:); 

%-----------initial total H2O, CH4 and salt masses-------------------------

init_w=sum(vb.*phi.*dnw.*sw.*(1-cl)+vb.*phi.*dnh.*sh*Mh2o);     % initial total water mass M/L^3
init_g=sum(vb.*phi.*dnh.*sh.*Mch4) +sum(vb.*phi.*dnw.*sw.*(1-cl).*cm) + sum(vb.*phi.*dng.*sg);      % initial total gas mass term
init_c=sum(vb.*phi.*dnw.*sw.*cl);       % initial total hydrate mass
init_w=sum(init_w);
init_g=sum(init_g);
init_c=sum(init_c);

