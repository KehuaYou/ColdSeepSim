function Rt=calc_Rt(i,j)
global COL
global pw pg  pcgw dpth vb
global twx twx1  tgx tgx1  ttx
global dnw dng  sh_0  T T_0 sg_0
global dt
global INDC1
global Cp_g Cp_w  qt  dy dz  L
global qg qw
global Heat_Capacity_ini

T0=0;
L=50000*(1/(0.016+5.75*0.018));
Mh=0.1195;
Rt=0;
pg=pw+pcgw;

adv_opt=1;
adv_opt_g=1;

if i==1
    if adv_opt==1
        if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
            Rt=Rt + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*Cp_w*(T0+T(i+1,j));
            Rt=Rt - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*Cp_w*(T0+T(i+1,j));
        else
            Rt=Rt + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
            Rt=Rt - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
        end
    end
    if adv_opt_g==1
        if tgx(i,j)*(pg(i+1,j)-pg(i,j)) > tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))
            Rt=Rt + tgx(i,j)*(pg(i+1,j)-pg(i,j))*dng(i+1,j)*Cp_g*(T0+T(i+1,j));
            Rt=Rt - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i+1,j)*Cp_g*(T0+T(i+1,j));
        else
            Rt=Rt + tgx(i,j)*(pg(i+1,j)-pg(i,j))*dng(i,j)*Cp_g*(T0+T(i,j));
            Rt=Rt - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i,j)*Cp_g*(T0+T(i,j));
        end
    end


    Rt=Rt + ttx(i,j)*(T(i+1,j)-T(i,j));




elseif i==COL
    if adv_opt==1
        if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
            Rt=Rt - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
            Rt=Rt + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
        else
            Rt=Rt - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*Cp_w*(T0+T(i-1,j));
            Rt=Rt + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*Cp_w*(T0+T(i-1,j));
        end
    end
    if adv_opt_g==1
        if tgx(i-1,j)*(pg(i,j)-pg(i-1,j)) > tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
            Rt=Rt - tgx(i-1,j)*(pg(i,j)-pg(i-1,j))*dng(i,j)*Cp_g*(T0+T(i,j));
            Rt=Rt + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i,j)*Cp_g*(T0+T(i,j));
        else
            Rt=Rt - tgx(i-1,j)*(pg(i,j)-pg(i-1,j))*dng(i-1,j)*Cp_g*(T0+T(i-1,j));
            Rt=Rt + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i-1,j)*Cp_g*(T0+T(i-1,j));
        end
    end

    Rt=Rt - ttx(i-1,j)*(T(i,j)-T(i-1,j));




else
    if adv_opt==1
        if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
            Rt=Rt + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*Cp_w*(T0+T(i+1,j));
            Rt=Rt - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*Cp_w*(T0+T(i+1,j));
        else
            Rt=Rt + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
            Rt=Rt - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
        end
    end
    if adv_opt_g==1
        if tgx(i,j)*(pg(i+1,j)-pg(i,j)) > tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))
            Rt=Rt + tgx(i,j)*(pg(i+1,j)-pg(i,j))*dng(i+1,j)*Cp_g*(T0+T(i+1,j));
            Rt=Rt - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i+1,j)*Cp_g*(T0+T(i+1,j));
        else
            Rt=Rt + tgx(i,j)*(pg(i+1,j)-pg(i,j))*dng(i,j)*Cp_g*(T0+T(i,j));
            Rt=Rt - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i,j)*Cp_g*(T0+T(i,j));
        end
    end

    Rt=Rt + ttx(i,j)*(T(i+1,j)-T(i,j));

    if adv_opt==1
        if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
            Rt=Rt - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
            Rt=Rt + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*Cp_w*(T0+T(i,j));
        else
            Rt=Rt - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*Cp_w*(T0+T(i-1,j));
            Rt=Rt + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*Cp_w*(T0+T(i-1,j));
        end
    end
    if adv_opt_g==1
        if tgx(i-1,j)*(pg(i,j)-pg(i-1,j)) > tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
            Rt=Rt - tgx(i-1,j)*(pg(i,j)-pg(i-1,j))*dng(i,j)*Cp_g*(T0+T(i,j));
            Rt=Rt + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i,j)*Cp_g*(T0+T(i,j));
        else
            Rt=Rt - tgx(i-1,j)*(pg(i,j)-pg(i-1,j))*dng(i-1,j)*Cp_g*(T0+T(i-1,j));
            Rt=Rt + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i-1,j)*Cp_g*(T0+T(i-1,j));
        end
    end

    Rt=Rt - ttx(i-1,j)*(T(i,j)-T(i-1,j));

end


Rt=-Rt + vb(i,j)*Heat_Capacity_ini(i,j)*(T(i,j)-T_0(i,j))/dt;


if (INDC1(i,j)==0 && sh_0(i,j)>0) || (INDC1(i,j)==1 && sg_0(i,j)>0)
    rh=Reaction_rate2(i,j);
    Rt=Rt-vb(i,j)*Mh*rh*L;
end


if i==COL
    Rt=Rt-dy(i,j)*dz(i,j)*qt-dy(i,j)*dz(i,j)*qg/(86400*365)*Cp_g*(T0+T(i,j))-dy(i,j)*dz(i,j)*qw/(86400*365)*Cp_w*(T0+T(i,j));
end
end