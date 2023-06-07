% Parameters
% k_in = 0.0197841057774864;%5.5e-5*360;
% k_olig_inc = 0.00538780048253155;%0.000015*360;
% k_olig_sep = 5.15448688362887e-06;%1.4e-8*360;
% k_clear_Abeta = 0.0198269961322812;%5.5e-5*360;
% k_onPP = 0.001*360;
% k_off_ma0 = 2290*k_onPP;
% k_off_ma1 = 67.3*k_onPP;
% k_plaque_inc = 2.52574188641094e-05;%7e-8*360;
% k_plaque_sep = 7.43970462839549e-09;%7e-11*360;
% k_clear_olig = 8.22167060462444e-06;%2.2e-8*360;
% k_clear_P = 1.11763490095714e-06;%4.41e-9*360;
% k_onPD = 0.001*360;
% k_off_ma2 = 1.79*k_onPD;
% k_synth_FcR = 0.0681207502332796;%(0.0000503/0.261)*360;
% k_clear_FcR = 0.0699452562807786;%0.000193*360;
% k_onPF = 0.001*360;
% k_offPF = 0.12*k_onPF;
% k_ADCP = 1.2904651505731;%0.0036*360;
% clearance = 0.00396;%1.1e-5*360;
% k_mAb_transport_back = 1.152;%0.0032*360;
% k_mAb_transport = 0.000576;%1.6e-6*360;
% k_mAbcomplex_clear = 5.4e-05;%1.5e-7*360;

k_in = 0.0197989883065441;%5.5e-5*360;
k_olig_inc = 0.00539992099177097;%0.000015*360;
k_olig_sep = 4.9994835431463e-06;%1.4e-8*360;
k_clear_Abeta = 0.0198002865326525;%5.5e-5*360;
k_onPP = 0.001*360;
k_off_ma0 = 2290*k_onPP;
k_off_ma1 = 67.3*k_onPP;
k_plaque_inc = 2.53440125750486e-05;%7e-8*360;
k_plaque_sep = 2.47973974755386e-08;%7e-11*360;
k_clear_olig = 7.66302511635033e-06;%2.2e-8*360;
k_clear_P = 1.51164568674074e-06 ;%4.41e-9*360;
k_onPD = 0.001*360;
k_off_ma2 = 1.79*k_onPD;
k_synth_FcR = 0.0693788627098226;%(0.0000503/0.261)*360;
k_clear_FcR = 0.069454611873771;%0.000193*360;
k_onPF = 0.001*360;
k_offPF = 0.12*k_onPF;
k_ADCP = 1.29601272340594;%0.0036*360;
clearance = 1.1e-5*360;
k_mAb_transport_back = 0.0032*360;
k_mAb_transport = 1.6e-6*360;
k_mAbcomplex_clear = 1.5e-7*360;

% Time range
tspan = 0:24:(364*24*1.5);

% Dose list
dose_list = 0:(24*14):(24*364*1.5);

% Initial conditions
initAbeta = 2;
initOlig = 4;
initPlaque = 39;
initFcR = 1;
initmAb = 0;
initAbetamAb = 0;
initOligmAb = 0;
initPlaquemAb = 0;
initOligmAbFcR = 0;
initPlaquemAbFcR = 0;
initPlasmamAb = 0;

initial_conditions = [initAbeta; initOlig; initPlaque; initFcR; initPlasmamAb; initmAb; initAbetamAb; initOligmAb; initPlaquemAb; initOligmAbFcR; initPlaquemAbFcR];
options = odeset('MaxStep',0.1);
[t,y] = ode15s(@(t,y) ODEs(t, y, k_in, k_olig_inc, k_olig_sep, k_clear_Abeta, ...
    k_onPP, k_off_ma0, k_off_ma1, k_plaque_inc, k_plaque_sep, k_clear_olig, ...
    k_clear_P, k_onPD, k_off_ma2, k_synth_FcR, k_clear_FcR, k_onPF, k_offPF, ...
    k_ADCP, clearance, k_mAb_transport_back, k_mAb_transport, k_mAbcomplex_clear, ...
    dose_list),tspan,initial_conditions, options);

%plot(t,y(:,5),'-o')
plaque_level = y(:,3)+y(:,11)+y(:,9);
percentage_mid = (((39)-plaque_level(53*7))/(39))*100
percentage_end = (((39)-plaque_level(end))/(39))*100
%plot(t,(y(:,3)+y(:,11)+y(:,9)),'-o')