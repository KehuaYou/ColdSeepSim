function Rc=calc_Rc(i,j)
global COL
global pw  dpth vb
global twx twx1
global sw sw_0 dnw dnw_0 cl cl_0
global phi phi_0 dt
global qw   dz  INDC2  dy
global Dsaltx  Cl_Bottom

Rc=0;
if i==1
    tmp=twx(i,j)*(pw(i+1,j)-pw(i,j)) - twx1(i,j)*(dpth(i+1,j)-dpth(i,j));
    if tmp > 0
        Rc=Rc + tmp*dnw(i+1,j)*cl(i+1,j);
    else
        Rc=Rc + tmp*dnw(i,j)*cl(i,j);
    end
    Rc=Rc + Dsaltx(i,j)*2*phi(i+1,j)*sw(i+1,j)*phi(i,j)*sw(i,j)/(phi(i+1,j)*sw(i+1,j)+phi(i,j)*sw(i,j))*(dnw(i+1,j)+dnw(i,j))/2*(cl(i+1,j)-cl(i,j));

elseif i==COL
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rc=Rc - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*cl(i,j);
        Rc=Rc + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*cl(i,j);
    else
        Rc=Rc - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*cl(i-1,j);
        Rc=Rc + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*cl(i-1,j);
    end
    Rc=Rc - Dsaltx(i-1,j)*2*phi(i,j)*sw(i,j)*phi(i-1,j)*sw(i-1,j)/(phi(i,j)*sw(i,j)+phi(i-1,j)*sw(i-1,j))*(dnw(i,j)+dnw(i-1,j))/2*(cl(i,j)-cl(i-1,j));

else
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rc=Rc + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*cl(i+1,j);
        Rc=Rc - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*cl(i+1,j);
    else
        Rc=Rc + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*cl(i,j);
        Rc=Rc - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*cl(i,j);
    end
    Rc=Rc + Dsaltx(i,j)*2*phi(i+1,j)*sw(i+1,j)*phi(i,j)*sw(i,j)/(phi(i+1,j)*sw(i+1,j)+phi(i,j)*sw(i,j))*(dnw(i+1,j)+dnw(i,j))/2*(cl(i+1,j)-cl(i,j));
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rc=Rc - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*cl(i,j);
        Rc=Rc + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*cl(i,j);
    else
        Rc=Rc - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*cl(i-1,j);
        Rc=Rc + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*cl(i-1,j);
    end
    Rc=Rc - Dsaltx(i-1,j)*2*phi(i,j)*sw(i,j)*phi(i-1,j)*sw(i-1,j)/(phi(i,j)*sw(i,j)+phi(i-1,j)*sw(i-1,j))*(dnw(i,j)+dnw(i-1,j))/2*(cl(i,j)-cl(i-1,j));
end


if INDC2(i,j)==1
    Rc=-Rc + ...
        vb(i,j)*phi(i,j)*dnw(i,j)*cl(i,j)*(1-sw_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*cl(i,j)*(dnw(i,j)-dnw_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*dnw_0(i,j)*(cl(i,j)-cl_0(i,j))/dt + ...
        vb(i,j)*sw_0(i,j)*dnw_0(i,j)*cl_0(i,j)*(phi(i,j)-phi_0(i,j))/dt;

else
    Rc=-Rc + ...
        vb(i,j)*phi(i,j)*dnw(i,j)*cl(i,j)*(sw(i,j)-sw_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*cl(i,j)*(dnw(i,j)-dnw_0(i,j))/dt + ...
        vb(i,j)*phi(i,j)*sw_0(i,j)*dnw_0(i,j)*(cl(i,j)-cl_0(i,j))/dt + ...
        vb(i,j)*sw_0(i,j)*dnw_0(i,j)*cl_0(i,j)*(phi(i,j)-phi_0(i,j))/dt;
end


if i==COL
    Rc=Rc-dy(i,j)*dz(i,j)*qw*Cl_Bottom/(1-Cl_Bottom)/(365*24*3600);
end

