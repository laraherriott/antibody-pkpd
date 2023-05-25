% Parameters
k_in = 5.5e-5;
k_olig_inc = 0.000015;
k_olig_sep = 1.4e-8;
k_clear_Abeta = 5.5e-5;
k_onPP = 0.001;
k_off_ma0 = 2290*k_onPP;
k_off_ma1 = 67.3*k_onPP;
k_plaque_inc = 7e-8;
k_plaque_sep = 7e-11;
k_clear_olig = 2.2e-8;
k_clear_P = 4.41e-9;
k_onPD = 0.001;
k_off_ma2 = 1.79*k_onPD;
k_synth_FcR = 0.0000503/0.261;
k_clear_FcR = 0.000193;
k_onPF = 0.001;
k_offPF = 0.12*k_onPF;
k_ADCP = 0.0036;
clearance = 0.00003;
k_mAb_transport_back = 0.0032;
k_mAb_transport = 1.6e-6;
k_mAbcomplex_clear = 1.5e-7;

% Time range
tspan = [0 (100*24*360)];

% Dose list
dose_list = 0:(24*14*360):(24*100*360);

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

[t,y] = ode45(@(t,y) ODEs(t, y, k_in, k_olig_inc, k_olig_sep, k_clear_Abeta, ...
    k_onPP, k_off_ma0, k_off_ma1, k_plaque_inc, k_plaque_sep, k_clear_olig, ...
    k_clear_P, k_onPD, k_off_ma2, k_synth_FcR, k_clear_FcR, k_onPF, k_offPF, ...
    k_ADCP, clearance, k_mAb_transport_back, k_mAb_transport, k_mAbcomplex_clear, ...
    dose_list),tspan,initial_conditions);
plot(t,y(:,5),'-o')
plot(t,(y(:,3)+y(:,11)+y(:,9)),'-o')