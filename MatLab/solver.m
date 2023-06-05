% Parameters
k_in =  5.5e-5*360;
k_olig_inc = 0.000015*360;
k_olig_sep = 1.4e-8*360;
k_clear_Abeta = 5.5e-5*360;
k_onPP = 0.001*360;
k_off_ma0 = 2290*k_onPP*360;
k_off_ma1 = 67.3*k_onPP*360;
k_plaque_inc = 7e-8*360;
k_plaque_sep = 7e-11*360;
k_clear_olig = 2.2e-8*360;
k_clear_P = 4.41e-9*360;
k_onPD = 0.001*360;
k_off_ma2 = 1.79*k_onPD*360;
k_synth_FcR = 0.0000503/0.261*360;
k_clear_FcR = 0.000193*360;
k_onPF = 0.001*360;
k_offPF = 0.12*k_onPF*360;
k_ADCP = 0.0036*360; %149.999985269331;
% clearance = 1.05e-5*360;
% k_mAb_transport_back = 0.0032*360;
% k_mAb_transport = 1.6e-6*360;
% k_mAbcomplex_clear = 1.5e-7*360;
clearance = 9.89185809639386e-05;%1.1e-5*360;
k_mAb_transport_back = 0.00310679015912649;%0.0032*360;
k_mAb_transport = 5.85746995234228e-06;%1.6e-6*360;
k_mAbcomplex_clear = 9.49706353503532e-07;%1.5e-7*360;

% Time range
tspan = 0:24:(364*24*1.5);

% Dose list
dose_list = 0:(24*14):(24*364*1.5);

% Initial conditions
initAbeta = 0.2;
initOlig = 370;
initPlaque = 5500;
initFcR = 1;
initmAb = 0;
initAbetamAb = 0;
initOligmAb = 0;
initPlaquemAb = 0;
initOligmAbFcR = 0;
initPlaquemAbFcR = 0;
initPlasmamAb = 0;

initial_conditions = [initAbeta; initOlig; initPlaque; initFcR; initPlasmamAb; initmAb; initAbetamAb; initOligmAb; initPlaquemAb; initOligmAbFcR; initPlaquemAbFcR];
options = odeset('MaxStep',0.01);
[t,y] = ode15s(@(t,y) ODEs(t, y, k_in, k_olig_inc, k_olig_sep, k_clear_Abeta, ...
    k_onPP, k_off_ma0, k_off_ma1, k_plaque_inc, k_plaque_sep, k_clear_olig, ...
    k_clear_P, k_onPD, k_off_ma2, k_synth_FcR, k_clear_FcR, k_onPF, k_offPF, ...
    k_ADCP, clearance, k_mAb_transport_back, k_mAb_transport, k_mAbcomplex_clear, ...
    dose_list),tspan,initial_conditions, options);

%plot(t,y(:,5),'-o')
plaque_level = y(:,3)+y(:,11)+y(:,9);
percentage_mid = (((5500)-plaque_level(53*7))/(5500))*100
percentage_end = (((5500)-plaque_level(end))/(5500))*100
%plot(t,(y(:,3)+y(:,11)+y(:,9)),'-o')