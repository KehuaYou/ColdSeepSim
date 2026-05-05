function Rh=calc_Rh(i,j)
global  sg_0
global phi phi_0 dt sh sh_0 dnh
global vb  INDC1

Mh=0.1195;  % molecular weight of methane hydrate

Rh=0;
Rh=Rh+vb(i,j)*phi(i,j)*dnh*(sh(i,j)-sh_0(i,j))/dt + ...
    vb(i,j)*sh_0(i,j)*dnh*(phi(i,j)-phi_0(i,j))/dt;

if (INDC1(i,j)==0 && sh_0(i,j)>0) || (INDC1(i,j)==1 && sg_0(i,j)>0)
    rh=Reaction_rate2(i,j);
    Rh=Rh-vb(i,j)*Mh*rh;
end



