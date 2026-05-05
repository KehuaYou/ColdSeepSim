clear
clear global
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global COL ROW
global INDC1 INDC2  INDC2_old
global cnstx cnsty  vb phi    kx ky
global pw pw_0 pcgw pcgw_0
global phi_0
global sw sw_0 sh sh_0 sg sg_0
global cl cl_0 cm cm_0
global dnw dnw_0 dng dng_0
global vsw vsw_0 vsg vsg_0
global dnh
global krw krg
global T T_0
global Mh2o Mch4 Msalt
global p_salinity T_salinity  salinity
global lambda
global dt t_flag  kx_old
global total_stress k_temp
global dx dy dz
global Fracture w0 C_frac S_frac Pcr  kx_ini  pcgw_ini
global twx twx1 tgx tgx1

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      Load  Initializatoin File             %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Option-1: run from initialization file
Initialization;         % load initial file
string = ('t0yr');      % save the initial results
save(string);           % save the initial results
timestep=1;

% Option-2: run from saved .mat file
% load('t8920yr.mat')      % load initial input from .mat file
% timestep=timestep+1;     % the start timestep=current timestep+1
% string_transient=('t8920yr.mat');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      time discretizaton            %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt_ini=86400*365*2;         % largest or initial time step
dt_cut_inner=0.5;           % when cut time step, dt=dt*dt_cut
dt_cut_outter=0.5;          % when we go back to previous saved timestep, dt_ini=dt_ini*dt_cut_outter
dt_minimum=1;               % the smallest timestep is 0.1 year
N=200000;                   % total number of timestep we will run
N_save=10;                  % we save our results every N_save timestep
N_show=50;                  % we show our results every N_save timestep

if timestep==1
    time=zeros(N+1,1);          % initializing time for each step
    time(1)=0;                  % we start at time=zero
    timestep=1;
    flux_w = zeros(N,1);        % seafloor water volumetric flux in m3/m2/sec
    flux_w_mass =zeros(N,1);    % seafloor water mass flux in kg/m2/sec
    flux_g = zeros(N,1);        % seafloor gas volumetric flux in m3/m2/sec
    flux_g_mass = zeros(N,1);   % seafloor gas mass flux in kg/m2/sec
    total_g=zeros(N,1);         % total amount of gas within the system, kg
    total_g(1)=init_g;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Main Subroutine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s

max_iter=5; % Maximum number of new_raphson iterations
go_back=0;  % flags to move on to next time step
break_flag=0;   % flags to kill the simulation

while timestep<=N
    disp(timestep)
    iteration=0;
    phase_check=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%       Save the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pw_0=pw;
    T_0=T;
    cm_0=cm;
    cl_0=cl;
    sw_0=sw;
    sh_0=sh;
    sg_0=sg;
    pcgw_0=pcgw;
    phi_0=phi;
    dnw_0=dnw;
    vsw_0=vsw;
    dng_0=dng;
    vsg_0=vsg;
    krw_old=krw;
    krg_old=krg;
    kx_old=kx;
    ky_old=ky;
    INDC1_old=INDC1;
    INDC2_old=INDC2;
    lambda_old=lambda;
    cnstx_old=cnstx;
    cnsty_old=cnsty;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % at new time step, we return to original dt
    dt=dt_ini;
    time(timestep+1)=time(timestep)+dt;
    cycle=0;
    loop=1;
    while loop==1
        go_back=0;
        iteration=iteration+1;
        newton_raphson_Kinetic;

        nan_flag=0; 
        for i=1:COL
            for j=1:ROW
                if isnan(pw(i,j)) || isnan(T(i,j))   % if the Jacobian matrix blows up, giving NAN values
                    nan_flag=1;
                    phase_check=1;
                    iteration=max_iter+1;   % get out of the loop
                    t_flag=1;               % decrease the timestep
                end
            end
        end
        if nan_flag==0
            ppt_update;
        end
        if iteration>max_iter 
            phase_check=1;
        end

        if phase_check==1
            if nan_flag==0
                cycle=cycle+1;
                %------------------- Primary Variable Switch --------------
                for i=1:COL
                    for j=1:ROW
                        if INDC2(i,j) == 2
                            INDC2(i,j)=3;
                            go_back=go_back+0;
                        end
                        if INDC2(i,j) == 1
                            if cm(i,j) > methane_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016
                                INDC2(i,j)=2;
                                go_back=go_back+1;
                            end
                        end

                        if INDC2(i,j)==3 && (sg(i,j)<-1e-4 || sh(i,j)<-1e-4)
                            INDC2(i,j)=2;
                            go_back=go_back+1;
                        end

                        if INDC2(i,j)==3 && sg(i,j)<-1e-4 && sh(i,j)<-1e-4
                            INDC2(i,j)=1;
                            go_back=go_back+1;
                        end
                    end
                end
                %-------------------------------------------------------
                if cycle>5
                    t_flag=1;
                end
                cm=max(0,cm);
            end

            if t_flag==1
                dt=dt_cut_inner*dt;
                disp(dt/86400/365)
                time(timestep+1)=time(timestep)+dt;
                INDC1=INDC1_old;
                INDC2=INDC2_old;
                t_flag=0;
                go_back=go_back+1;
                cycle=0;
            end
            if go_back>0
                pw=pw_0;
                cm=cm_0;
                cl=cl_0;
                sw=sw_0;
                sh=sh_0;
                sg=sg_0;
                pcgw=pcgw_0;
                phi=phi_0;
                dnw=dnw_0;
                vsw=vsw_0;
                dng=dng_0;
                vsg=vsg_0;
                T=T_0;

                krw=krw_old;
                krg=krg_old;
                kx=kx_old;
                ky=ky_old;
                lambda=lambda_old;
                cnstx=cnstx_old;
                cnsty=cnsty_old;

                loop=1;
                iteration=0;
                phase_check=0;
            else
                loop=0;
            end % end of if go_back==1
        end  % end of phase_check=1
    end % end of while loop==1

    if loop==0
        %%%%%%%%%%%%%%%%%%%%%%% Update phase zone %%%%%%%%%%%%%%%%%%%%%%%
        for i=1:COL
            for j=1:ROW
                if cl(i,j) < interpolation2(p_salinity,T_salinity,salinity, pw(i,j)/1e6,T(i,j))
                    INDC1(i,j)=1;
                else
                    INDC1(i,j)=0;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%% Seafloor fluid flux %%%%%%%%%%%%%%%%%%%%%%%
        % -------- Water flux
        flux_w(timestep,1)= twx(1,1)*(pw(2,1)-pw(1,1)) - twx1(1,1)*(dpth(2,1)-dpth(1,1));
        flux_w_mass(timestep,1)= flux_w(timestep,1)*dnw(2,1)*(1-cl(2,1));

        % --------- Gas flux
        flux_g(timestep,1) = tgx(1,1)*(pw(2,1)+pcgw(2,1)-pw(1,1)-pcgw(1,1)) - tgx1(1,1)*(dpth(2,1)-dpth(1,1));
        flux_g_mass(timestep,1) = flux_g(timestep,1)*dng(2,1);

        % --------- Total mass of methane within the system
        total_g(timestep)=sum(vb.*phi.*dnh.*sh.*Mch4) +sum(vb.*phi.*dnw.*sw.*(1-cl).*cm) + sum(vb.*phi.*dng.*sg);

        %%%%%%%%%%%%%%%%%% Fracture propogation %%%%%%%%%%%%%%%%%%%%%%%%%%
        eff_stress=total_stress-pw;
        if pw(COL,1)> total_stress(COL,1)+Pcr
            Fracture(COL,1)=1;
        end
        for i=1:COL-1
            if pw(i,1)> total_stress(i,1)+Pcr || pw(i+1,1)>total_stress(i,1)+Pcr
                Fracture(i,1)=1;
            end
        end
        Fracture(1,1)=Fracture(2,1);

        % --------- Permeability and capillary pressure update
        w_frac = w0*exp(-C_frac.*max(eff_stress,0));
        w_frac(1)=w0;
        k_temp=(1-sh).^2.*w_frac.^3./(6*S_frac).*Fracture+kx_ini.*(1-Fracture);
        pcgw=(2*0.074./(w_frac.*sw)).*Fracture+pcgw_ini.*(1-Fracture);
        kx=k_temp;

        % --------- Flow operator update
        for i=1:COL
            for j=1:ROW
                if i<COL
                    cnstx(i,j)=kx(i+1,j)*(dy(i,j)*dz(i,j)/(dx(i,j)*0.5+dx(i+1,j)*0.5));
                else
                    cnstx(i,j)=kx(i,j)*(dy(i,j)*dz(i,j)/(dx(i,j)*0.5+dx(i-1,j)*0.5));
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%% Show and Save outputs %%%%%%%%%%%%%%%%%%%%%%%
        
        % ---------- Show outputs

        if ~mod(timestep,N_show)

            % Outputs: depth below seafloor, temperature, total stress,
            % water pressure, effective stresss, fracture distribution,
            % fracture size, permeability, water saturation, gas
            % saturation, hydrate saturation, salinity, dissolved methane
            % concentration, index for phase zone, index for phase number

            junk=[  (dpth(:,1)-seafloor)./1e3 T(:,1)  total_stress(:,1)./1e6 ...
                pw(:,1)./1e6 eff_stress(:,1)./1e6 Fracture(:, 1) ...
                w_frac(:,1).*1e6 kx(:,1).*1e12  sw(:,1)  sg(:,1)  sh(:,1) cl(:,1) ...
                cm(:,1)  INDC1(:,1)  INDC2(:,1) ]
        end

        % ---------- Save ouputs
        
        if ~mod(timestep,N_save)
            string = ['t',num2str(time(timestep+1)/(86400*365)),'yr','.mat'];
            string_transient=string;
            save(string)
        end
        timestep=timestep+1;
    end

    if dt<dt_minimum
        if timestep<=N_save
            load('t0ky.mat');
            timestep=1;
        else
            load(string_transient);
        end
        disp('we have gone back to previous saved timestep')
        dt_ini=dt_ini*dt_cut_outter;

        if timestep>1
            timestep=timestep+1;
        end
    end
end

final_w=sum(sum(vb.*phi.*dnw.*sw.*(1-cl)+vb.*phi.*dnh.*sh*Mh2o));       % neglect the mass of methane in water when calculate water mass
final_g=sum(sum(vb.*phi.*dnh.*sh.*Mch4)) + sum(sum(vb.*phi.*dng.*sg)) + sum(sum(vb.*phi.*dnw.*sw.*(1-cl).*cm));
final_c=sum(sum(vb.*phi.*dnw.*sw.*cl));
final_w=sum(final_w);
final_g=sum(final_g);
final_c=sum(final_c);














