global COL ROW
global dnw
global dng
global krw krg
global pw
global sg sh
global T
global cl cm
global Msalt

for i=1:COL
    for j=1:ROW
        if sg(i,j)>0
            cm(i,j)=methane_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016;
        elseif sh(i,j)>0
            cm(i,j)=hydrate_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016;
        end

        dng(i,j)=gas_density(pw(i,j)/1e6,T(i,j));

        dnw(i,j)=brine_density(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)),cm(i,j)/0.016);

        [krw(i,j), krg(i,j)] = calc_relative_perm(i,j);
    end
end

