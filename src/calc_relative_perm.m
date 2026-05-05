function [krw_calc, krg_calc] = calc_relative_perm(i,j)
global sw sh
global wn gn
global swr sgr
global INDC2

if sw(i,j)>swr
    krw_calc=((sw(i,j)-swr)/(1-swr-sgr))^wn;
else
    krw_calc=0;
end

if INDC2(i,j) > 1
    sg_temp = 1 -sw(i,j)-sh(i,j);
    if sg_temp>= sgr
        krg_calc=((sg_temp-sgr)/(1-swr-sgr))^gn;
    else
        krg_calc=0;
    end
else
    krg_calc=0;
end

end