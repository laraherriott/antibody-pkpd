function err=funLSQ(p, y_observed, k_onPP, k_onPD, k_onPF, k_offPF, k_off_ma0, k_off_ma1, k_off_ma2)
k_in = p(1);
k_olig_inc = p(2);
k_olig_sep = p(3);
k_clear_Abeta = p(4);
k_plaque_inc = p(5);
k_plaque_sep = p(6);
k_clear_olig = p(7);
k_clear_P = p(8);
k_synth_FcR = p(9);
k_clear_FcR = p(10);
k_ADCP = p(11);
clearance = p(12);
k_mAb_transport_back = p(13);
k_mAb_transport = p(14);
k_mAbcomplex_clear = p(15);

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

initial_conditions = [initAbeta; initOlig; initPlaque; initFcR; initmAb; initAbetamAb; initOligmAb; initPlaquemAb; initOligmAbFcR; initPlaquemAbFcR; initPlasmamAb];
options=odeset('MaxStep',0.1);
[t,y] = ode15s(@(t,y) ODEs(t, y, k_in, k_olig_inc, k_olig_sep, k_clear_Abeta, ...
    k_onPP, k_off_ma0, k_off_ma1, k_plaque_inc, k_plaque_sep, k_clear_olig, ...
    k_clear_P, k_onPD, k_off_ma2, k_synth_FcR, k_clear_FcR, k_onPF, k_offPF, ...
    k_ADCP, clearance, k_mAb_transport_back, k_mAb_transport, k_mAbcomplex_clear, ...
    dose_list),tspan,initial_conditions, options);

plaque = y(:,3)+y(:,9)+y(:,11);
plaque_initial = plaque(1);
plaque_fall = zeros(length(plaque)-1,1);

for i = (1:1:length(plaque)-1)
    plaque_fall(i,1) = ((plaque_initial-plaque(i))/plaque_initial)*100;
end

plaque_change = zeros(length(plaque)-1,1);
plaque_change(53*7,1) = plaque_fall(53*7,1);
plaque_change(end,1) = plaque_fall(end,1);
plaque_error = plaque_change - y_observed(:,1);

ab = zeros(length(plaque)-1,1);

for i = (1:1:length(plaque)-1)
    if i <= 100
        ab(i,1) = y(i,5);
    end
end
ab_error = ab - y_observed(:,2);
ab_error2 = max(ab) - max(y_observed(:,2));

% Time range
tspan = 0:24:(24*365);

% Dose list
dose_list = 0:0:0;

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

initial_conditions = [initAbeta; initOlig; initPlaque; initFcR; initmAb; initAbetamAb; initOligmAb; initPlaquemAb; initOligmAbFcR; initPlaquemAbFcR; initPlasmamAb];
options=odeset('MaxStep',0.1);
[t,y2] = ode15s(@(t,y) ODEs(t, y, k_in, k_olig_inc, k_olig_sep, k_clear_Abeta, ...
    k_onPP, k_off_ma0, k_off_ma1, k_plaque_inc, k_plaque_sep, k_clear_olig, ...
    k_clear_P, k_onPD, k_off_ma2, k_synth_FcR, k_clear_FcR, k_onPF, k_offPF, ...
    k_ADCP, clearance, k_mAb_transport_back, k_mAb_transport, k_mAbcomplex_clear, ...
    dose_list),tspan,initial_conditions, options);

plaque2 = y2(:,3)+y2(:,9)+y2(:,11);
plaque_initial2 = plaque2(1);
plaque_fall2 = zeros(length(plaque2)-1,1);

for i = (1:1:length(plaque2)-1)
    plaque_fall2(i,1) = ((plaque_initial2-plaque2(i))/plaque_initial2)*100;
end

plaque_change2 = zeros(length(plaque2)-1,1);
plaque_change2(end,1) = plaque_fall2(end,1);
plaque_change2(numel(plaque_error)) = 0;
plaque_error2 = plaque_change2;

% ab_error = y(:,5) - y_observed(1:31,2);
% ab_error(numel(plaque_error)) = 0;
err=ab_error+10*plaque_error+100*(plaque_error2);
end

