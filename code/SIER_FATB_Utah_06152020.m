
clear all
close all

%% LOOP

%% Read in data
datamatrix=[];
datamatrix = importdata('../data/state44.csv');
    %(1) ccases
    %(2) hospitalized
    %(3) death
    %(4) aa
    %(5) lockdown
    %(6) a0
    %(7) a1
    %(8) a2
    %(9) a3
    %(10) lambda
    %(11) pop


%% Step 1: Data preparation
    %Basic estimation period
    data.N=mean(datamatrix(:,11)); %Population
    TT=length(datamatrix(:,4)) - sum(isnan(datamatrix(:,4))); %Time horizon
    Tend = length(datamatrix(:,4));

    %calibrated or pre-estimated
    data.sigma = 1/5.2;  %rate at which exposed become symptomatic
    data.gamma = 1/5;    %1/5 rate at which infectious become "resolving"
    data.theta = 1/14;   %1/14 rate at which disease is resolved (transit to recovered or dead)

    %Lockdown policy calibration
    data.eta_L = 1/3;
    data.lambda_P=mean(datamatrix(:,10));

    data.Tstart_lock = find(datamatrix(:,5), 1, 'first');
    data.Tend_lock = find(datamatrix(:,5), 1, 'last');

    %Behavioral response parameters
        data.a0 = mean(datamatrix(:,6));      %constant
        data.a1 = mean(datamatrix(:,7));     %confirmed cases
        data.a2 = 0* mean(datamatrix(:,8));    %hospitalized
        data.a3 = mean(datamatrix(:,9));      %deaths

%         data.a0 = 1.2;      %constant
%         data.a1 = -0.05;     %confirmed cases
%         data.a2 = 0;    %hospitalized
%         data.a3 = 0.01;      %deaths

%         data.a0 = 1.177;      %constant
%         data.a1 = -0.0857;     %confirmed cases
%         data.a2 = 0;    %hospitalized
%         data.a3 = 0.086;      %deaths
%         data.lambda_P=-0.093;

%% Loop for out of sample estimation

Nlearners = 28;

kk = -(Nlearners-1):0;

for l = 1:Nlearners;

data.T = TT +kk(l);

    %load data
        data.O_m = datamatrix(1:data.T,1);
        data.H_m = datamatrix(1:data.T,2);
        data.F_m = datamatrix(1:data.T,3);
        data.a_m = datamatrix(1:data.T,4);




%% Step 2: Estimation

    %Set up estimation
    lb(1,1) = 6;         ub(1,1) = 6;         %R_0 virality
    lb(2,1) = 0;         ub(2,1) = 1;         %tauR: random testing rate
    lb(3,1) = 0;         ub(3,1) = 0;        %tauS1: symptom-based testing of infectious people
    lb(4,1) = 0;         ub(4,1) = 1;        %tauS2: symptom-based testing of resolving people
    lb(5,1) = 0;         ub(5,1) = 1000000;    %initial exposed (latent)
    lb(6,1) = 0.01;         ub(6,1) = 0.8;        %alpha: probability of not displaying symptoms
    lb(7,1) = 0;         ub(7,1) = 1;        %delta: death rate
    lb(8,1) = 1;         ub(8,1) = 10;        %psi: infection curvature on activities
    lb(9,1) = 0;         ub(9,1) = 1;        %gamma: transition infectious-> resolving
    lb(10,1) =1/12;         ub(10,1) = 1/12;       %theta: transition resolving->end
    lb(11,1) =1/5;         ub(11,1) = 1/5;       %sigma: transition exposed -> infectious; from: https://www.cdc.gov/coronavirus/2019-ncov/hcp/clinical-guidance-management-patients.html#:~:text=The%20incubation%20period%20for%20COVID,CoV%2D2%20infection.

    lb(12,1) =0;         ub(12,1) = 1000000;     %Initial I^US
    lb(13,1) =0;         ub(13,1) = 1000000;     %Initial I^UA
    lb(14,1) =0;         ub(14,1) = 1000000;     %Initial R^US
    lb(15,1) =0;         ub(15,1) = 1000000;     %Initial R^UA





    %initialize estimation
    theta0=[];
    theta0(1)=6;
    theta0(2)=0;
    theta0(3)=0;
    theta0(4)=0.001;
    theta0(5)=100;
    theta0(6)=0.8;
    theta0(7)=0.08;
    theta0(8)=2;
    theta0(9)=1/5;
    theta0(10)=1/10;
    theta0(11)=1/5;

    theta0(12)=1;
    theta0(13)=1;
    theta0(14)=1;
    theta0(15)=1;


    %define function objective handle
    inSEE_SEIR_FATB_fh = @(theta0)inSEE_SEIR_FATB12(theta0,data);

        %inSEE_BSEIRFAT_fh([theta_true(1), theta_true(2),theta_true(3), theta_true(4)])
    sol_zknitro=[];
    [sol_zknitro, fval_zknitro, exitflag,output,lambda] = knitro_nlp(inSEE_SEIR_FATB_fh, theta0, [], [], [], [], lb, ub, [], [], []);
    %sol_zknitro(8)=2;
%% Step 3: save key estimation results

sol_zknitro_mat(l,:) = sol_zknitro;
exitflags(l,1) = exitflag;



%% Step: set up prediction file
%
% para.sigma = sol_zknitro(11);
% para.gamma = sol_zknitro(9);
% para.theta = sol_zknitro(10);
% para.delta = sol_zknitro(7);
% para.psi = sol_zknitro(8)
% para.R_0 = sol_zknitro(1);
% hp = sol_zknitro(6);
% para.N=data.N;
% para.T=data.T;
% para.Tfor= length(datamatrix(:,1))- data.T;
% para.Tpolicy=30;
% para.tauR= sol_zknitro(2);
% para.tauS1= sol_zknitro(3);
% para.tauS2= sol_zknitro(4);
% para.lambda_P=data.lambda_P;
% para.Tstart_lock = data.Tstart_lock;
% para.Tend_lock= data.Tend_lock;
% para.eta_L= data.eta_L;
% para.a0 = data.a0;
% para.a1 = data.a1;
% para.a2 = data.a2;
% para.a3 = data.a3;
% para.R_US_1 = sol_zknitro(14);
% para.R_UA_1 = sol_zknitro(15);
% para.I_US_1 = sol_zknitro(12);
% para.I_UA_1 = sol_zknitro(13);
% para.E_1= sol_zknitro(5);
% para.alpha = sol_zknitro(6);
%
% predmat = SIER_FATB_predictions1(para);






%% Step 4: display estimation results (ans save more results)

sigma = sol_zknitro(11);  %rate at which exposed become symptomatic
gamma = sol_zknitro(9);    %rate at which infectious become "resolving"
theta = sol_zknitro(10);   %rate at which disease is resolved (transit to recovered or dead)
delta = sol_zknitro(7);   %death rate 0.08=0.8% death rate
R_0 = sol_zknitro(1);        %initial virality
beta0=R_0*gamma;         % initial beta
alpha = sol_zknitro(6);    %Probability of NOT displaying symptoms

N=data.N; %Population
T=data.T; %Time horizon

Tfor = length(datamatrix(:,1))- data.T;
Tpolicy=30;
%% Policy parameters

    %Policy 1: Testing & quarantining
        tauR = sol_zknitro(2);      %random testing rate
        tauS1 = sol_zknitro(3);    %testing of symptomatic infectious people
        tauS2 = sol_zknitro(4);    %testing of symptomatic infectious people

        %Tstart_test = 60;
        %kappat = kappa.*[zeros(Tstart_test,1); ones(T-Tstart_test,1)];

    %Policy 2: Forced lockdown
        lambda_P=data.lambda_P; %forced lockdown, e.g. lambda=0.5 locks down 40% of the workforce
        Tstart_lock = data.Tstart_lock;
        Tend_lock = data.Tend_lock;
        eta_L = data.eta_L;

   %Policy 3: Information policies
        a0 = data.a0;      %constant
        a1 = data.a1;     %confirmed cases
        a2 = data.a2;    %hospitalized
        a3 = data.a3;      %deaths


%% Initialize model

S = zeros(T,1);
R_DS = zeros(T,1);
R_DA = zeros(T,1);
R_D = zeros(T,1);
R_U = zeros(T,1);
R_A = zeros(T,1);
R_S = zeros(T,1);
R_US = zeros(T,1);
R_UA = zeros(T,1);
I_DS = zeros(T,1);
I_DA = zeros(T,1);
I_D = zeros(T,1);
I_US = zeros(T,1);
I_UA = zeros(T,1);
I_U = zeros(T,1);
I_S = zeros(T,1);
I_U = zeros(T,1);
E_D = zeros(T,1);
E = zeros(T,1);
O = zeros(T,1);
F = zeros(T,1);
C = zeros(T,1);
C_D = zeros(T,1);
a = zeros(T,1);
lambdat=zeros(T-1+Tfor+Tpolicy,1);

ls = S;
betat = S;
    betat(1) = beta0;
a(1) = a0;

%% Initialization of disease
R_DS(1) = 1;                          %one initial case (hospitalized)
R_DA(1) = 0;                          %one initial case
    R_D(1) = R_DS(1)+R_DA(1);

R_US(1) = sol_zknitro(14);    %implied non-hospitalized resolving cases
R_UA(1) = sol_zknitro(15);
    R_U(1) = R_US(1)+R_UA(1);

    R_S(1) =R_DS(1)+R_US(1);
    R_A(1) = R_UA(1)+R_DA(1);

I_DS(1) = 0;
I_DA(1) = 0;
    I_D(1) = I_DS(1)+I_DA(1);

I_US(1) = sol_zknitro(12);
I_UA(1) = sol_zknitro(13);
    I_U(1) = I_US(1)+I_UA(1);

    I_S(1) = I_DS(1)+I_US(1);
    I_A(1) = I_DA(1)+I_UA(1);

E_D(1) = 0;
%E(1) = (1/sigma)*(I_DS(1)+I_DA(1)+I_US(1)+I_UA(1));
E(1) = sol_zknitro(5);

O(1) = R_D(1);

S(1) = N - E(1) - I_D(1) - I_U(1) - R_D(1) - R_U(1);

F(1) = 0;
C(1) = 0;

%% Calculate time path

for t = 1:T-1+Tfor+Tpolicy;
    %% Public health module
        %Susceptiables
        S(t+1) = S(t) -betat(t)*S(t)*(I_U(t))/N;

        %Exposed
        E(t+1) = E(t) + betat(t)*S(t)*(I_U(t))/N - sigma*E(t) - tauR*E(t);
        E_D(t+1) = E_D(t) + tauR*E(t) - sigma*E_D(t);

        %Infected - undetected
        I_UA(t+1) = I_UA(t) + alpha*sigma*E(t) -gamma*I_UA(t) - tauR*I_UA(t);       %asymptomatic
        I_US(t+1) = I_US(t) + (1-alpha)*sigma*E(t) -gamma*I_US(t) - tauS1*I_US(t) - tauR*I_US(t);   %symptomatic
        I_U(t+1) = I_UA(t+1) + I_US(t+1);

        %Infected - detected
        I_DA(t+1) = I_DA(t) + alpha*sigma*E_D(t) - gamma*I_DA(t) + tauR*I_UA(t);        %asymptomatic
        I_DS(t+1) = I_DS(t) + (1-alpha)*sigma*E_D(t) - gamma*I_DS(t) + tauS1*I_US(t)+tauR*I_US(t);    %symptomatic
        I_D(t+1) = I_DA(t+1) + I_DS(t+1);

        %Resolving - undetected
        R_UA(t+1) = R_UA(t) + gamma*I_UA(t) - theta*R_UA(t) - tauR*R_UA(t);   %asymptomatic
        R_US(t+1) = R_US(t) + gamma*I_US(t) - theta*R_US(t) - tauS2*R_US(t)- tauR*R_US(t);   %symptomatic
            R_U(t+1) = R_UA(t+1) + R_US(t+1);

        %Resolving - detected
        R_DA(t+1) = R_DA(t) + gamma*I_DA(t) - theta*R_DA(t) + tauR*R_UA(t);   %asymptomatic
        R_DS(t+1) = R_DS(t) + gamma*I_DS(t) - theta*R_DS(t) + tauS2*R_US(t)+tauR*R_US(t);   %symptomatic

            R_D(t+1) = R_DA(t+1) + R_DS(t+1);

            R_S(t+1) = R_US(t+1) + R_DS(t+1);
            R_A(t+1) = R_UA(t+1) + R_DA(t+1);


        %Fatalities
        F(t+1) = F(t) + delta*(theta*R_US(t)+theta*R_DS(t));

        %ReCovered
        C(t+1) = C(t) + theta*R_UA(t) + (1-delta)*theta*R_US(t);
        C_D(t+1) = C_D(t) + theta*R_DA(t) + (1-delta)*theta*R_DS(t);

        %Hospitalizations
        H(t+1) = R_DS(t+1);

        %Observable cases
        O(t+1) = E_D(t+1) + I_D(t+1) + H(t+1);




        %% Behavioral Response module

        %Lockdown policy effect

    if Tstart_lock==[];
        lambdat(t) = 0;
    else
            if t < Tstart_lock
            lambdat(t) = 0;
            elseif t >= Tstart_lock & t < Tend_lock
            lambdat(t) = data.lambda_P*(1-exp(-eta_L*(t-Tstart_lock)));
            elseif t >= Tend_lock
            lambdat(t) = lambdat(Tend_lock-1)*exp(-eta_L*(t-Tend_lock));
            end
    end



        %Endogenous avoidance
        %a(t) = exp( log(a_0) - xi1*( 100*(O(t+1))/N )^eta1 - xi2*( 100*R_DS(t)/N )^eta2 - xi3*( 100*(F(t+1)-F(t))/N )^eta3 -lambdat(t) );

        %a(t+1) = a(t) - lambdat(t) - 0.1*a(t);


        %grO(t+1) = 2*(O(t+1)-O(t))/(O(t+1)+O(t));
        %grH(t+1) = 2*(H(t+1)-H(t))/(H(t+1)+H(t));
        %grF(t+1) = 2*(F(t+1)-F(t))/(F(t+1)+F(t));
        %a(t+1) = a0 + a1*log(1+O(t+1)) + a2*grO(t+1) + a3*log(1+H(t+1)) + a4*grH(t+1) + a5*log(1+F(t+1)) + a6*grF(t+1);

        a(t+1) = a0 + a1*log(1+O(t+1))+ a2*log(1+H(t+1)) +  a3*log(1+F(t+1)) + lambdat(t);

        %% Economic module

        %Employment share effects due to sickness
        ls_sick(t) = ( S(t)+E(t)+ I_UA(t) + R_UA(t) +R_DA(t) + C(t) )/(N-F(t));

        %Total employment as share of initial employment
        ls(t) = ls_sick(t)*a(t);

        % Avoidance health effect
        betat(t+1) = beta0*a(t)^sol_zknitro(8);

end

%% Save results for out of sample predictions
        O_mat(l,:) = O(1:(data.T+Tfor+Tpolicy));
        F_mat(l,:) = F(1:(data.T+Tfor+Tpolicy));
        a_mat(l,:) = a(1:(data.T+Tfor+Tpolicy));
end






%% Set up the ensemble

ensembledata.Tfor_min = TT;
ensembledata.Tfor_max = data.T+Tfor;
ensembledata.O_mat = O_mat(:,ensembledata.Tfor_min:ensembledata.Tfor_max);
ensembledata.F_mat = F_mat(:,ensembledata.Tfor_min:ensembledata.Tfor_max);
ensembledata.a_mat = a_mat(:,ensembledata.Tfor_min:ensembledata.Tfor_max);

ensembledata.Odat = datamatrix(ensembledata.Tfor_min:ensembledata.Tfor_max,1)';
ensembledata.Fdat = datamatrix(ensembledata.Tfor_min:ensembledata.Tfor_max,3)';
ensembledata.adat = datamatrix(ensembledata.Tfor_min:ensembledata.Tfor_max,4)';


wgt0 = ones(1,Nlearners)/Nlearners;

%% Define function handle
ensemble_cve_fh = @(wgt0)ensemble_cve(wgt0,ensembledata);

%% Define constraints

    lb =[];
    ub=[];

    lb = zeros(length(wgt0),1);
    ub = ones(length(wgt0),1);
    Aeq = ones(1,length(wgt0));
    beq=1;

%% solve ensemble problem
    [wgt_opt, fval_zknitro, exitflag,output,lambda] = knitro_nlp(ensemble_cve_fh, wgt0, [], [], Aeq, beq, lb, ub, [], [], []);

    ensemblepath_O = sum(kron(wgt_opt',ones(1,length(O_mat))).*O_mat);
    ensemblepath_F = sum(kron(wgt_opt',ones(1,length(O_mat))).*F_mat);
    ensemblepath_a = sum(kron(wgt_opt',ones(1,length(O_mat))).*a_mat);

%% Plot out results

 close all

figure(1)
hax=axes;
    h1 = plot(1:(data.T+Tfor), datamatrix(:,1), '*');
    hold on
    h2 = plot(1:(data.T+Tfor+Tpolicy), O_mat(l,:), '--');
    h3 = plot(1:(data.T+Tfor+Tpolicy), O_mat(2,:), '--');
    h4 = plot(1:(data.T+Tfor+Tpolicy), O_mat(3,:), '--');
    h5 = plot(1:(data.T+Tfor+Tpolicy), O_mat(4,:), '--');
    h6 = plot(1:(data.T+Tfor+Tpolicy), O_mat(5,:), '--');
    h7 = plot(1:(data.T+Tfor+Tpolicy), O_mat(6,:), '--');
    h8 = plot(1:(data.T+Tfor+Tpolicy), O_mat(7,:), '--');
    h9 = plot(1:(data.T+Tfor+Tpolicy), O_mat(8,:), '--');
    h10= plot(1:(data.T+Tfor+Tpolicy), O_mat(9,:), '--');
    h11 = plot(1:(data.T+Tfor+Tpolicy), O_mat(10,:), '--');
    h12 = plot(1:(data.T+Tfor+Tpolicy), mean(O_mat), '--');
    h13 = plot(1:(data.T+Tfor+Tpolicy), ensemblepath_O, '-');
    h14 = line([data.T data.T],get(hax,'YLim'), 'LineWidth',1,'Color',[1 0 0]);
    h15 = line([data.T+kk(1) data.T+kk(1)],get(hax,'YLim'), 'LineWidth',1,'Color',[0.5 0.5 0.5]);
    hold off
    fh = figure(1); % returns the handle to the figure object
    axis([0 115 0 15000])
        xlabel('Days since first documented infection','fontsize',25); ylabel('Observed cases','fontsize',25);
        grid on;
    colormap gray;
    %    lgd = legend({'Data','Training sample to \tau-12','Training sample to \tau-11','Training sample to \tau-10', 'Training sample to \tau-9', 'Training sample to \tau-8','Training sample to \tau-7','Training sample to \tau-6','Training sample to \tau-5','Training sample to \tau-4','Training sample to \tau-3','Training sample to \tau-2','Training sample to \tau-1','Training sample to \tau','Median Ensemlble','Optimal Ensemble', 'End of estimation sample \tau'}, 'Location', 'NorthWest','NumColumns',2);
    %lgd.FontSize = 5;
    set(fh, 'color', 'white'); % sets the color to white
    set(gca, 'Box', 'off', 'fontsize', 16);
    set(h1, 'Marker', '*', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 5.0);
    set(h12, 'LineWidth', 2, 'Color', [1 0 0]);
    set(h13, 'LineWidth', 2, 'Color', [0 0 1]);


figure(2)
hax=axes;
    h1 = plot(1:(data.T+Tfor), datamatrix(:,3), '*');
    hold on
    h2 = plot(1:(data.T+Tfor+Tpolicy), F_mat(l,:), '--');
    h3 = plot(1:(data.T+Tfor+Tpolicy), F_mat(2,:), '--');
    h4 = plot(1:(data.T+Tfor+Tpolicy), F_mat(3,:), '--');
    h5 = plot(1:(data.T+Tfor+Tpolicy), F_mat(4,:), '--');
    h6 = plot(1:(data.T+Tfor+Tpolicy), F_mat(5,:), '--');
    h7 = plot(1:(data.T+Tfor+Tpolicy), F_mat(6,:), '--');
    h8 = plot(1:(data.T+Tfor+Tpolicy), F_mat(7,:), '--');
    h9 = plot(1:(data.T+Tfor+Tpolicy), F_mat(8,:), '--');
    h10= plot(1:(data.T+Tfor+Tpolicy), F_mat(9,:), '--');
    h11 = plot(1:(data.T+Tfor+Tpolicy), F_mat(10,:), '--');
    h12 = plot(1:(data.T+Tfor+Tpolicy), median(F_mat), '--');
    h13 = plot(1:(data.T+Tfor+Tpolicy), ensemblepath_F, '-');
    h14 = line([data.T data.T],get(hax,'YLim'), 'LineWidth',1,'Color',[1 0 0]);
    h15 = line([data.T+kk(1) data.T+kk(1)],get(hax,'YLim'), 'LineWidth',1,'Color',[0.5 0.5 0.5]);
    hold off
    fh = figure(2); % returns the handle to the figure object
        xlabel('Days since first documented infection','fontsize',25); ylabel('Cumulative fatalities','fontsize',25);
        grid on;
    colormap gray;
        axis([0 115 0 200])
        %lgd = legend({'Data','Training sample to \tau-10', 'Training sample to \tau-9', 'Training sample to \tau-8','Training sample to \tau-7','Training sample to \tau-6','Training sample to \tau-5','Training sample to \tau-4','Training sample to \tau-3','Training sample to \tau-2','Training sample to \tau-1','Training sample to \tau','End of estimation sample \tau'}, 'Location', 'NorthWest','NumColumns',2);
    %lgd.FontSize = 5;
    set(fh, 'color', 'white'); % sets the color to white
    set(gca, 'Box', 'off', 'fontsize', 16);
    set(h1, 'Marker', '*', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 5.0);
    set(h12, 'LineWidth', 2, 'Color', [1 0 0]);
    set(h13, 'LineWidth', 2, 'Color', [0 0 1]);



figure(3)
hax=axes;
    h1 = plot(1:(data.T+Tfor), datamatrix(:,4), '*');
    hold on
    h2 = plot(1:(data.T+Tfor+Tpolicy), a_mat(l,:), '--');
    h3 = plot(1:(data.T+Tfor+Tpolicy), a_mat(2,:), '--');
    h4 = plot(1:(data.T+Tfor+Tpolicy), a_mat(3,:), '--');
    h5 = plot(1:(data.T+Tfor+Tpolicy), a_mat(4,:), '--');
    h6 = plot(1:(data.T+Tfor+Tpolicy), a_mat(5,:), '--');
    h7 = plot(1:(data.T+Tfor+Tpolicy), a_mat(6,:), '--');
    h8 = plot(1:(data.T+Tfor+Tpolicy), a_mat(7,:), '--');
    h9 = plot(1:(data.T+Tfor+Tpolicy), a_mat(8,:), '--');
    h10= plot(1:(data.T+Tfor+Tpolicy), a_mat(9,:), '--');
    h11 = plot(1:(data.T+Tfor+Tpolicy), a_mat(10,:), '--');
    h12 = plot(1:(data.T+Tfor+Tpolicy), median(a_mat), '--');
    h13 = plot(1:(data.T+Tfor+Tpolicy), ensemblepath_a, '-');
    h14 = line([data.T data.T],get(hax,'YLim'), 'LineWidth',1,'Color',[1 0 0]);
    h15 = line([data.T+kk(1) data.T+kk(1)],get(hax,'YLim'), 'LineWidth',1,'Color',[0.5 0.5 0.5]);
    hold off
    fh = figure(3); % returns the handle to the figure object
        xlabel('Days since first documented infection','fontsize',25); ylabel('Mobility relative to 2019','fontsize',25);
        grid on;
    colormap gray;
    axis([0 115 0.5 1.2])
        %lgd = legend({'Data','Training sample to \tau-10', 'Training sample to \tau-9', 'Training sample to \tau-8','Training sample to \tau-7','Training sample to \tau-6','Training sample to \tau-5','Training sample to \tau-4','Training sample to \tau-3','Training sample to \tau-2','Training sample to \tau-1','Training sample to \tau','End of estimation sample \tau'}, 'Location', 'NorthWest','NumColumns',2);
    %lgd.FontSize = 5;
    set(fh, 'color', 'white'); % sets the color to white
    set(gca, 'Box', 'off', 'fontsize', 16);
    set(h1, 'Marker', '*', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 5.0);
    set(h12, 'LineWidth', 2, 'Color', [1 0 0]);
    set(h13, 'LineWidth', 2, 'Color', [0 0 1]);

%
%      print(figure(1),'-dpdf', 'fit_ccases')
%      print(figure(2),'-dpdf', 'fit_fatalities')
%      print(figure(3),'-dpdf', 'fit_mobility')
%
%
    %
% %     figure(2)
% %     hax=axes;
% %     h1 = plot(1:(data.T+Tfor), datamatrix(:,2), '*');
% %     hold on
% %     h2 = plot(1:(data.T+Tfor), H(1:(data.T+Tfor)), '-');
% %     h3 = line([data.T data.T],get(hax,'YLim'), 'LineWidth',1,'Color',[1 0 0]);
% %     hold off
% %     fh = figure(2); % returns the handle to the figure object
% %         xlabel('Days since first documented infection','fontsize',25); ylabel('Hospitalized','fontsize',25);
% %         grid on;
% %     colormap gray;
% %     legend('Data', 'Estimated hospital bed demand','End of estimation sample', 'Location', 'NorthWest')
% %     set(fh, 'color', 'white'); % sets the color to white
% %     set(gca, 'Box', 'off', 'fontsize', 16);
% %     set(h1, 'Marker', '*', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 5.0);
% %     set(h2, 'LineWidth', 2, 'Color', [0 0 1]);
% %     title([{sprintf('%s',statename{l})} {'Hospital bed demand'}]);
%
%     figure(2)
%     hax=axes;
%     h1 = plot(1:(data.T+Tfor), datamatrix(:,3), '*');
%     hold on
%     h2 = plot(1:(data.T+Tfor), F(1:(data.T+Tfor)), '-');
%     h3 = line([data.T data.T],get(hax,'YLim'), 'LineWidth',1,'Color',[1 0 0]);
%     hold off
%     fh = figure(2); % returns the handle to the figure object
%         xlabel('Days since first documented infection','fontsize',25); ylabel('Cumulative deaths','fontsize',25);
%         grid on;
%     colormap gray;
%     legend('Data', 'Estimated fatality path', 'Location', 'NorthWest')
%     set(fh, 'color', 'white'); % sets the color to white
%     set(gca, 'Box', 'off', 'fontsize', 16);
%     set(h1, 'Marker', '*', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 5.0);
%     set(h2, 'LineWidth', 2, 'Color', [0 0 1]);
%     title([{sprintf('%s',statename{l})} {'Fatalities'}]);
%
% figure(3)
%     h1 = plot(1:data.T, datamatrix(1:data.T,4), '*');
%     hold on
%     h2 = plot(1:data.T, a(1:data.T), '-');
%     hold off
%     fh = figure(3); % returns the handle to the figure object
%         xlabel('Days since first documented infection','fontsize',25); ylabel('Activities (relative to 2019)','fontsize',25);
%         grid on;
%     colormap gray;
%     legend('Data', 'Estimated activity path', 'Location', 'NorthWest')
%     set(fh, 'color', 'white'); % sets the color to white
%     set(gca, 'Box', 'off', 'fontsize', 16);
%     set(h1, 'Marker', '*', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 5.0);
%     set(h2, 'LineWidth', 2, 'Color', [0 0 1]);
%     title([{sprintf('%s',statename{l})} {'Economic Activities'}]);
%
%         figure(5)
%     h1 = plot(datamatrix(1:data.T,3)/datamatrix(1:data.T,2));
%     hold on
%     h2 = plot(F(1:(data.T+Tfor))/H(1:(data.T+Tfor))');
%     hold off

%     %save graphs
%     cd 'E:\Nowcasting COVID\SIR\MJ\results\graphs'
%
%     print(figure(1),'-dpdf', fit_fig1{l})
%     print(figure(2),'-dpdf', fit_fig2{l})
%     print(figure(3),'-dpdf', fit_fig3{l})1
%     print(figure(4),'-dpdf', fit_fig4{l})
%
%     %Saving more results
%     Results(l,1).O = O;
%     Results(l,1).H = H;
%     Results(l,1).F = F;
%     Results(l,1).a = a;
%     Results(l,1).betat = length(datamatrix(:,1))- data.T;
%     Results(l,1).Tpolicy = Tpolicy;
%
%     Results(l,1).TR_A = betat.*S.*I_UA./data.N;
%     Results(l,1).TR_S = betat.*S.*I_US./data.N;



%% Extract estimates
%
% cd 'E:\Nowcasting COVID\SIR\MJ\Paper1'
%
%     %Calibrated parameters
%     Tcalibrated = table( ['R_0  ';'theta';'sigma'], [ sol_zknitro_mat(10,1) ; sol_zknitro_mat(10,10) ; sol_zknitro_mat(10,11) ] );
%     writetable(Tcalibrated)
%
%     %Estimates
%     Testimates = table(['tauR ';'tauS ';'alpha'; 'delta'; 'psi  '; 'gamma'], [ sol_zknitro_mat(10,2) ; sol_zknitro_mat(10,4) ; sol_zknitro_mat(10,6) ; sol_zknitro_mat(10,7) ; sol_zknitro_mat(10,8) ; sol_zknitro_mat(10,9) ]);
%     writetable(Testimates)
%
%
% %% Calculate MSE
%
% for l = 1:10
%     inR2_O(l,1) = 1- sum( (O(1:TT+kk(l)) - datamatrix(1:TT+kk(l),1)).^2 )/sum( (datamatrix(1:TT+kk(l),1)).^2 );
%     inAVG_O(l,1) = sum(O(1:TT+kk(l)))/(TT+kk(l));
%     outR2_O(l,1) = 1- sum( (O(TT+kk(l):Tend) - datamatrix(TT+kk(l):Tend,1)).^2 )/sum( (datamatrix(TT+kk(l):Tend,1)).^2 );
%     outAVG_O(l,1) = sum(O(TT+kk(l):Tend))/(Tend-(TT+kk(l)));
%
%     inR2_F(l,1) = 1- sum( (F(1:TT+kk(l)) - datamatrix(1:TT+kk(l),3)).^2 )/sum( (datamatrix(1:TT+kk(l),3)).^2 );
%     inAVG_F(l,1) = sum(F(1:TT+kk(l)))/(TT+kk(l));
%     outR2_F(l,1) = 1- sum( (F(TT+kk(l):Tend) - datamatrix(TT+kk(l):Tend,3)).^2 )/sum( (datamatrix(TT+kk(l):Tend,3)).^2 );
%     outAVG_F(l,1) = sum(F(TT+kk(l):Tend))/(Tend-(TT+kk(l)));
%
%     inR2_a(l,1) = 1- sum( (a(1:TT+kk(l)) - datamatrix(1:TT+kk(l),4)).^2 )/sum( (datamatrix(1:TT+kk(l),4)).^2 );
%     inAVG_a(l,1) = sum(a(1:TT+kk(l)))/(TT+kk(l));
%     outR2_a(l,1) = 1- sum( (a(TT+kk(l):TT) - datamatrix(TT+kk(l):TT,4)).^2 )/sum( (datamatrix(TT+kk(l):TT,4)).^2 );
%     outAVG_a(l,1) = sum(a(TT+kk(l):Tend))/(Tend-(TT+kk(l)));
% end
%
%
% %plot(kk(1:10), out2_O)
% %plot(kk(1:10), inR2_O)
%
% Tmodelfit = table( ['inAVG_O '; 'outAVG_O'; 'inR2_O  '; 'outR2_O ' ; 'inAVG_F '; 'outAVG_F'; 'inR2_F  '; 'outR2_F ' ; 'inAVG_a '; 'outAVG_a'; 'inR2_a  '; 'outR2_a '], ...
%    [inAVG_O(10,1); outAVG_O(10,1); inR2_O(10,1); outR2_O(10,1); inAVG_F(10,1); outAVG_F(10,1); inR2_F(10,1); outR2_F(10,1); inAVG_a(10,1); outAVG_a(10,1); inR2_a(10,1); outR2_a(10,1)] );
% writetable(Tmodelfit)
%

% kk = [2,4,5,7,8,9,12,13,14,15];
%
% for k = 1:10
%     coefvar(k,1) = std(sol_zknitro_mat(:,kk(k)))/mean(sol_zknitro_mat(:,kk(k)))
% end
