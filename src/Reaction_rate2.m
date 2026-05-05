function Rate=Reaction_rate2(i,j)
global  dt dnh  sh sg phi
global phi_0 sg_0 dng_0 sh_0
global pw T cl Msalt dnw
global p_salinity T_salinity salinity
global INDC1
global S_frac
global hydrate_on

if hydrate_on == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                                                                  %%%%%%
    %%%                         constant parameters                      %%%%%%
    %%%                                                                  %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mm=0.016;           % molar weight of methane
    Mh=0.1195;          % molar weight of hydrate
    Dm=1e-16;           % methane diffusion coefficient through hydrate skin

    k_kinetic=0.5875e-11*1e-3;      

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                                                                  %%%%%%
    %%%      methane hydrate/gas solubility in liquid water              %%%%%%
    %%%                                                                  %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C_gw=methane_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016;
    C_hw=hydrate_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016;
    C_gw=dnw(i,j)*(1-cl(i,j))*C_gw/Mm;
    C_hw=dnw(i,j)*(1-cl(i,j))*C_hw/Mm;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                                                                  %%%%%%
    %%% formation rate limited by methane diffusion through hydrate skin %%%%%%
    %%% dissocation rate limited by Kim-Bishinoi kinetic model           %%%%%%
    %%%                                                                  %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sh_temp=sh(i,j);
    if sh_temp<1e-4 && sg(i,j)>0
        sh_temp=1e-4;
    end

    if INDC1(i,j)==1 % hydrate formation

        if sh(i,j)<0.9999
            Dm=1e-16;
        else
            Dm=1e-20;
        end
        A_gw=4*(sg(i,j)+sh(i,j))/S_frac;
        rate_reaction=A_gw^2/(sh_temp*phi(i,j))*Dm*(C_gw-C_hw);

        if isnan(rate_reaction)
            rate_reaction=0;
        end

    else  % hydrate dissociation

        As=4*(sg(i,j)+sh(i,j))/S_frac;
        fugacity_equi=1e6*interpolation_pressure(p_salinity,T_salinity,salinity,cl(i,j),T(i,j));
        fugacity_current=pw(i,j);
        rate_reaction=k_kinetic*As*(fugacity_current-fugacity_equi);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                                                                  %%%%%%
    %%%      rate limited by the amount of local material                %%%%%%
    %%%                                                                  %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if INDC1(i,j)==1
        rate_material=phi_0(i,j)*sg_0(i,j)*dng_0(i,j)/Mm/dt;
    else
        rate_material=-phi_0(i,j)*sh_0(i,j)*dnh/Mh/dt;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                                                                  %%%%%%
    %%%    final = the smaller one between reation and material          %%%%%%
    %%%                                                                  %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if abs(rate_reaction)>abs(rate_material)
        Rate=rate_material;
    else
        Rate=rate_reaction;
    end

else

    Rate = 0;

end

