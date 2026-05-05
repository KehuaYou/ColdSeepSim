function Rg=calc_Rg(i,j)
global COL
global pw pcgw  dpth
global tgx tgx1  twx twx1
global dnw dnw_0 sw sw_0 dng dng_0  sg_0 sh sh_0 cl cl_0 cm cm_0
global phi phi_0 dt
global INDC1 INDC2 INDC2_old
global vb
global qg  dz
global dy
global Dmethanex

Mm=0.016;
Rg=0;
if i==1
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rg=Rg + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*(1-cl(i+1,j))*cm(i+1,j);
        Rg=Rg - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*(1-cl(i+1,j))*cm(i+1,j);
    else
        Rg=Rg + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
        Rg=Rg - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
    end
    if tgx(i,j)*(pw(i+1,j)+pcgw(i+1,j)-pw(i,j)-pcgw(i,j)) > tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rg=Rg + tgx(i,j)*(pw(i+1,j)-pw(i,j))*dng(i+1,j) + tgx(i,j)*(pcgw(i+1,j)-pcgw(i,j))*dng(i+1,j);
        Rg=Rg - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i+1,j);
    else
        Rg=Rg + tgx(i,j)*(pw(i+1,j)-pw(i,j))*dng(i,j) + tgx(i,j)*(pcgw(i+1,j)-pcgw(i,j))*dng(i,j);
        Rg=Rg - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i,j);
    end
    Rg=Rg + Dmethanex(i,j)*(phi(i+1,j)*sw(i+1,j)+phi(i,j)*sw(i,j))/2*(dnw(i+1,j)+dnw(i,j))/2*(cm(i+1,j)-cm(i,j));

elseif i==COL
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rg=Rg - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
        Rg=Rg + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
    else
        Rg=Rg - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j))*cm(i-1,j);
        Rg=Rg + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j))*cm(i-1,j);
    end
    if tgx(i-1,j)*(pw(i,j)+pcgw(i,j)-pw(i-1,j)-pcgw(i-1,j)) > tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rg=Rg - tgx(i-1,j)*(pw(i,j)-pw(i-1,j))*dng(i,j) - tgx(i-1,j)*(pcgw(i,j)-pcgw(i-1,j))*dng(i,j);
        Rg=Rg + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i,j);
    else
        Rg=Rg - tgx(i-1,j)*(pw(i,j)-pw(i-1,j))*dng(i-1,j) - tgx(i-1,j)*(pcgw(i,j)-pcgw(i-1,j))*dng(i-1,j);
        Rg=Rg + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i-1,j);
    end
    Rg=Rg - Dmethanex(i-1,j)*(phi(i,j)*sw(i,j)+phi(i-1,j)*sw(i-1,j))/2*(dnw(i,j)+dnw(i-1,j))/2*(cm(i,j)-cm(i-1,j));

else
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rg=Rg + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*(1-cl(i+1,j))*cm(i+1,j);
        Rg=Rg - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*(1-cl(i+1,j))*cm(i+1,j);
    else
        Rg=Rg + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
        Rg=Rg - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
    end
    if tgx(i,j)*(pw(i+1,j)+pcgw(i+1,j)-pw(i,j)-pcgw(i,j)) > tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rg=Rg + tgx(i,j)*(pw(i+1,j)-pw(i,j))*dng(i+1,j) + tgx(i,j)*(pcgw(i+1,j)-pcgw(i,j))*dng(i+1,j);
        Rg=Rg - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i+1,j);
    else
        Rg=Rg + tgx(i,j)*(pw(i+1,j)-pw(i,j))*dng(i,j) + tgx(i,j)*(pcgw(i+1,j)-pcgw(i,j))*dng(i,j);
        Rg=Rg - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i,j);
    end
    Rg=Rg + Dmethanex(i,j)*(phi(i+1,j)*sw(i+1,j)+phi(i,j)*sw(i,j))/2*(dnw(i+1,j)+dnw(i,j))/2*(cm(i+1,j)-cm(i,j));

    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rg=Rg - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
        Rg=Rg + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
    else
        Rg=Rg - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j))*cm(i-1,j);
        Rg=Rg + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j))*cm(i-1,j);
    end
    if tgx(i-1,j)*(pw(i,j)+pcgw(i,j)-pw(i-1,j)-pcgw(i-1,j)) > tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rg=Rg - tgx(i-1,j)*(pw(i,j)-pw(i-1,j))*dng(i,j) - tgx(i-1,j)*(pcgw(i,j)-pcgw(i-1,j))*dng(i,j);
        Rg=Rg + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i,j);
    else
        Rg=Rg - tgx(i-1,j)*(pw(i,j)-pw(i-1,j))*dng(i-1,j) - tgx(i-1,j)*(pcgw(i,j)-pcgw(i-1,j))*dng(i-1,j);
        Rg=Rg + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i-1,j);
    end
    Rg=Rg - Dmethanex(i-1,j)*(phi(i,j)*sw(i,j)+phi(i-1,j)*sw(i-1,j))/2*(dnw(i,j)+dnw(i-1,j))/2*(cm(i,j)-cm(i-1,j));


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%       change in mass accumulation term            %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if INDC2(i,j) == 1
    Rg=-Rg + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*(1-cl(i,j))*(cm(i,j)-cm_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*cm_0(i,j)*(cl_0(i,j)-cl(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(dnw(i,j)-dnw_0(i,j))/dt + ...
        vb(i,j)*dnw_0(i,j)*sw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(phi(i,j)-phi_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(1-sw_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*dng_0(i)*(0-sg_0(i,j))/dt;  % this line is added



elseif INDC2(i,j) == 2

    if INDC2_old(i,j)==1
        if INDC1(i,j)==1
            Rg=-Rg + ...
                vb(i,j)*phi(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)*(sw(i,j)-sw_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*(1-cl(i,j))*(cm(i,j)-cm_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*cm_0(i,j)*(cl_0(i,j)-cl(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(dnw(i,j)-dnw_0(i,j))/dt + ...
                vb(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(phi(i,j)-phi_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*dng(i,j)*(1-sw(i,j)-sh(i,j)-sg_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sg_0(i,j)*(dng(i,j)-dng_0(i,j))/dt + ...
                vb(i,j)*dng_0(i,j)*sg_0(i,j)*(phi(i,j)-phi_0(i,j))/dt;
        else
            Rg=-Rg + ...
                vb(i,j)*phi(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)*(sw(i,j)-sw_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*(1-cl(i,j))*(cm(i,j)-cm_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*cm_0(i,j)*(cl_0(i,j)-cl(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(dnw(i,j)-dnw_0(i,j))/dt + ...
                vb(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(phi(i,j)-phi_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*dng(i,j)*(1-sw(i,j)-sg_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sg_0(i,j)*(dng(i,j)-dng_0(i,j))/dt + ...
                vb(i,j)*dng_0(i,j)*sg_0(i,j)*(phi(i,j)-phi_0(i,j))/dt;
        end
    elseif INDC2_old(i,j)==3
        if INDC1(i,j) == 0   %%% sh=0
            Rg=-Rg + ...
                vb(i,j)*phi(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)*(sw(i,j)-sw_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*(1-cl(i,j))*(cm(i,j)-cm_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*cm_0(i,j)*(cl_0(i,j)-cl(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(dnw(i,j)-dnw_0(i,j))/dt + ...
                vb(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(phi(i,j)-phi_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*dng(i,j)*(1-sw(i,j)-sg_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sg_0(i,j)*(dng(i,j)-dng_0(i,j))/dt + ...
                vb(i,j)*dng_0(i,j)*sg_0(i,j)*(phi(i,j)-phi_0(i,j))/dt;
        else  %% sg=0
            Rg=-Rg + ...
                vb(i,j)*phi(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)*(sw(i,j)-sw_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*(1-cl(i,j))*(cm(i,j)-cm_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*cm_0(i,j)*(cl_0(i,j)-cl(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(dnw(i,j)-dnw_0(i,j))/dt + ...
                vb(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(phi(i,j)-phi_0(i,j))/dt+ ...
                vb(i,j)*phi(i,j)*dng(i,j)*(0-sg_0(i,j))/dt + ...
                vb(i,j)*phi(i,j)*sg_0(i,j)*(dng(i,j)-dng_0(i,j))/dt + ...
                vb(i,j)*dng_0(i,j)*sg_0(i,j)*(phi(i,j)-phi_0(i,j))/dt;

        end
    end

elseif INDC2(i,j) == 3

    Rg=-Rg + ...
        vb(i,j)*phi(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)*(sw(i,j)-sw_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*(1-cl(i,j))*(cm(i,j)-cm_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*dnw(i,j)*cm_0(i,j)*(cl_0(i,j)-cl(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(dnw(i,j)-dnw_0(i,j))/dt + ...
        vb(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)*(phi(i,j)-phi_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*dng(i,j)*(1-sw(i,j)-sh(i,j)-sg_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sg_0(i,j)*(dng(i,j)-dng_0(i,j))/dt + ...
        vb(i,j)*dng_0(i,j)*sg_0(i,j)*(phi(i,j)-phi_0(i,j))/dt;

end

if (INDC1(i,j)==0 && sh_0(i,j)>0) || (INDC1(i,j)==1 && sg_0(i,j)>0)
    rm=Reaction_rate2(i,j);
    Rg=Rg+vb(i,j)*Mm*rm;
end

%++++++++++++++++++++inject gas from below+++++++++++++++++++++++++
if i==COL
    Rg=Rg-dy(i,j)*dz(i,j)*qg/(365*24*3600);
end



