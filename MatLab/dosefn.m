function sol = dosefn(dose_list, t)
    infusion = (((((10*70)/1000/3.22)/147181.62))*1e9);
    f = 0.5;
    delta = 0.1;
  
    sol = 0;
    for n = dose_list
        if t >=(n-0.5) && t <= (n+(1.5))
            sol = ((infusion/2)/atan(1/delta))*(atan(sin(2*pi*(t-n-0.5)*f)/delta)) + (infusion/2);
        end
    end


end