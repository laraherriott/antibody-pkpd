% Parameters
k_in = 5.5e-5*360;
k_olig_inc = 0.000015*360;
k_olig_sep = 1.4e-8*360;
k_clear_Abeta = 5.5e-5*360;
k_onPP = 0.001*360;
k_off_ma0 = 2290*k_onPP;
k_off_ma1 = 67.3*k_onPP;
k_plaque_inc = 7e-8*360;
k_plaque_sep = 7e-11*360;
k_clear_olig = 2.2e-8*360;
k_clear_P = 4.41e-9*360;
k_onPD = 0.001*360;
k_off_ma2 = 1.79*k_onPD;
k_synth_FcR = (0.0000503/0.261)*360;
k_clear_FcR = 0.000193*360;
k_onPF = 0.001*360;
k_offPF = 0.12*k_onPF;
k_ADCP = 0.0036*360;
clearance = 1.1e-5*360;
k_mAb_transport_back = 0.0032*360;
k_mAb_transport = 1.6e-6*360;
k_mAbcomplex_clear = 1.5e-7*360;

% k_in =  0.0198012632058048;%5.5e-5*360;
% k_olig_inc = 0.00540954435109234;%0.000015*360;
% k_olig_sep = 1.28592552012923e-05;%1.4e-8*360;
% k_clear_Abeta = 0.019808619723581;%5.5e-5*360;
% k_onPP = 0.001*360;
% k_off_ma0 = 2290*k_onPP*360;
% k_off_ma1 = 67.3*k_onPP*360;
% k_plaque_inc = 2.51998542285385e-05;%7e-8*360;
% k_plaque_sep = 2.31705493570231e-06;%7e-11*360;
% k_clear_olig = 1.09606924358198e-05;%2.2e-8*360;
% k_clear_P = 1.58758599989717e-06;%4.41e-9*360;
% k_onPD = 0.001*360;
% k_off_ma2 = 1.79*k_onPD*360;
% k_synth_FcR = 0.069389154851369 ;%0.0000503/0.261*360;
% k_clear_FcR = 0.0694823966820438;%0.000193*360;
% k_onPF = 0.001*360;
% k_offPF = 0.12*k_onPF*360;
% k_ADCP = 128;%1.29599813749577;%0.0036*360;
% clearance = 0.00399982999733133; %0.00003*360;
% k_mAb_transport_back = 1.1520003314899 ; %0.0032*360;
% k_mAb_transport = 0.000576833653356959; %1.6e-6*360;
% k_mAbcomplex_clear = 5.45890853286003e-05; %1.5e-7*360;

% % Time range
% tspan = 0:360:(50*24*360);
% 
% % Dose list
% dose_list = 0:(24*14*360):(24*100*360);
% 
% % Initial conditions
% initAbeta = 0.2;
% initOlig = 370;
% initPlaque = 5500;
% initFcR = 1;
% initmAb = 0;
% initAbetamAb = 0;
% initOligmAb = 0;
% initPlaquemAb = 0;
% initOligmAbFcR = 0;
% initPlaquemAbFcR = 0;
% initPlasmamAb = 0;
% 
% initial_conditions = [initAbeta; initOlig; initPlaque; initFcR; initmAb; initAbetamAb; initOligmAb; initPlaquemAb; initOligmAbFcR; initPlaquemAbFcR; initPlasmamAb];
% 
% [t,y] = ode45(@(t,y) ODEs(t, y, k_in, k_olig_inc, k_olig_sep, k_clear_Abeta, ...
%     k_onPP, k_off_ma0, k_off_ma1, k_plaque_inc, k_plaque_sep, k_clear_olig, ...
%     k_clear_P, k_onPD, k_off_ma2, k_synth_FcR, k_clear_FcR, k_onPF, k_offPF, ...
%     k_ADCP, clearance, k_mAb_transport_back, k_mAb_transport, k_mAbcomplex_clear, ...
%     dose_list),tspan,initial_conditions);

% n_observations = length(t);
% scale = 0.1;
% noise = scale*randn(n_observations, 1);
% y_observed = y(:,3) + noise;


ab = csvread('../PK_fit_data.csv', 1, 1);
plaque_observed = zeros([364*1.5 1]);
plaque_observed(53*7, 1) = 67.5;
plaque_observed(end, 1) = 77.5; 
ab_observed = zeros([364*1.5 1]);
for i = (1:1:length(ab_observed))
    count=0;
    if i <= 100
        ab_observed(i, 1) = ab(i*24,1);
    end

end

% ab_decay = csvread('../../lec_decay_one_dose.csv',1,0);
% ab_d_obs = zeros([364*1.5 1]);
% count=1;
% for i = (1:1:length(ab_d_obs))
%     if ismember(i, ab_decay(:,1))
%         ab_d_obs(i, 1) = ab_decay(count,2);
%         count = count + 1;
%     end
% 
% end

y_observed = [plaque_observed, ab_observed];
%y_observed2 = ab_d_obs;

initial_estimate = [k_ADCP];
lb = [0];%, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
ub = [10];%[0.1, 0.01, 1e-5, 0.1, 1e-4, 1e-6, 1e-4, 1e-4, 0.1, 0.1, 10]; 

options = optimoptions('lsqnonlin', 'OptimalityTolerance',1e-8);

%[xlsqnonlin, errorlsqnonlin]  = lsqnonlin(@(p)funLSQ(p, y_observed,  k_in, k_olig_inc, k_olig_sep, k_clear_Abeta, k_plaque_inc, k_plaque_sep, k_clear_olig, k_clear_P, k_synth_FcR, k_clear_FcR, clearance, k_mAb_transport_back, k_mAb_transport, k_mAbcomplex_clear, k_onPP, k_onPD, k_onPF, k_offPF, k_off_ma0, k_off_ma1, k_off_ma2), initial_estimate, lb, ub, options)


problem = createOptimProblem('lsqnonlin', 'x0', initial_estimate, 'objective', @(p)funLSQ(p, y_observed, k_in, k_olig_inc, k_olig_sep, k_clear_Abeta, k_plaque_inc, k_plaque_sep, k_clear_olig, k_clear_P, k_synth_FcR, k_clear_FcR, clearance, k_mAb_transport_back, k_mAb_transport, k_mAbcomplex_clear, k_onPP, k_onPD, k_onPF, k_offPF, k_off_ma0, k_off_ma1, k_off_ma2), 'lb', lb, 'ub', ub);
%stpoints = RandomStartPointSet('NumStartPoints', 5, 'ArtificialBound', 100);

ms = MultiStart('Display', 'iter');
ms.TolFun = 1e-8;
[xmultinonlin,errormultinonlin, eflag, output] = run(ms, problem, 10)