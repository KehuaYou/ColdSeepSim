function Rw=calc_Rw(i,j)
global COL 
global pw  dpth
global twx twx1 
global sw sw_0 dnw dnw_0 cl cl_0
global phi phi_0  sh_0  sg_0
global  dt
global  INDC1 INDC2
global vb
global qw  dz  dy

Mw=0.018;
N_hydration=5.75;
Rw=0;
if i==1
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
    else
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*(1-cl(i,j));
    end


elseif i==COL
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*(1-cl(i,j));
    else
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
    end

else
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
    else
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*(1-cl(i,j));
    end
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*(1-cl(i,j));
    else
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%       change in mass accumulation term            %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if INDC2(i,j) == 1
    Rw=-Rw + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*(1-cl(i,j))*(dnw(i,j)-dnw_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*dnw_0(i,j)*(cl_0(i,j)-cl(i,j))/dt + ...
        vb(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*(phi(i,j)-phi_0(i,j))/dt+ ...
        vb(i,j)*phi(i,j)*dnw(i,j)*(1-cl(i,j))*(1-sw_0(i,j))/dt;

else

    Rw=-Rw + ...
        vb(i,j)*phi(i,j)*dnw(i,j)*(1-cl(i,j))*(sw(i,j)-sw_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*(1-cl(i,j))*(dnw(i,j)-dnw_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*dnw_0(i,j)*(cl_0(i,j)-cl(i,j))/dt + ...
        vb(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*(phi(i,j)-phi_0(i,j))/dt;
end

if (INDC1(i,j)==0 && sh_0(i,j)>0) || (INDC1(i,j)==1 && sg_0(i,j)>0)
    rw=N_hydration*Reaction_rate2(i,j);
    Rw=Rw+vb(i,j)*Mw*rw;
end

if i==COL
    Rw=Rw-dy(i,j)*dz(i,j)*qw/(365*24*3600);
end


