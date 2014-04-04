%% lkt_date3

function [dt]=lkt_dateindex(T,frac,start_y,start_m,index_s,index_e,sub_st)
    if start_m > frac; error('Invalid starting month: %d, with fraction: %d.',start_m,frac); end
    y = floor(start_y+(start_m-1)/frac+(0:frac*T-1)'./frac);
    f = mod((0:frac*T-1)'+start_m,frac); f2 = f == 0; f(f2) = frac;
    dt = [y f]; dt = dt(1:T,:); dt = dt(:,1) + (dt(:,2)-1)./frac;
end
