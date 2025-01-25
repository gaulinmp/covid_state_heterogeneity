%% Objective function for estimation





function output = inSEE_SEIR_FATB12(x, data)

%% Parameters



sigma = x(11);  %rate at which exposed become symptomatic
gamma = x(9);    %rate at which infectious become "resolving"
theta = x(10);   %rate at which disease is resolved (transit to recovered or dead)
delta = x(7);   %death rate 0.08=0.8% death rate
R_0 = x(1);        %initial virality
beta0=R_0*gamma;         % initial beta
alpha = x(6);    %Probability of NOT displaying symptoms

N=data.N; %Population
T=data.T; %Time horizon

%% Policy parameters

    %Policy 1: Testing & quarantining
        tauR = x(2);      %random testing rate
        tauS1 = x(3);    %testing of symptomatic infectious people
        tauS2 = x(4);    %testing of symptomatic infectious people        

        %Tstart_test = 60;
        %kappat = kappa.*[zeros(Tstart_test,1); ones(T-Tstart_test,1)];
        
    %Policy 2: Forced lockdown (calibrated)
        lambda_P=data.lambda_P; %forced lockdown, e.g. lambda=0.5 locks down 40% of the workforce    
        Tstart_lock = data.Tstart_lock;
        Tend_lock = data.Tend_lock;
        eta_L = data.eta_L;
        
        if Tstart_lock==[];
            Tstart_lock=0;
        end
       
   %Policy 3: Information policies
        a0 = data.a0;      %constant
        a1 = data.a1;     %confirmed cases
        a2 = data.a2;   %hospitalized
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
lambdat = zeros(T,1);

ls = S;
betat = S;
    betat(1) = beta0;
a(1) = a0;

%% Initialization of disease
R_DS(1) = 1;                          %one initial case (hospitalized)
R_DA(1) = 0;                          %one initial case
    R_D(1) = R_DS(1)+R_DA(1);

R_US(1) = x(14);    %implied non-hospitalized resolving cases
R_UA(1) = x(15);
    R_U(1) = R_US(1)+R_UA(1);

    R_S(1) =R_DS(1)+R_US(1);
    R_A(1) = R_UA(1)+R_DA(1);
    
I_DS(1) = 0;
I_DA(1) = 0;
    I_D(1) = I_DS(1)+I_DA(1);

I_US(1) = x(12);
I_UA(1) = x(13);
    I_U(1) = I_US(1)+I_UA(1);

    I_S(1) = I_DS(1)+I_US(1);
    I_A(1) = I_DA(1)+I_UA(1);
    
E_D(1) = 0;
%E(1) = (1/sigma)*(I_DS(1)+I_DA(1)+I_US(1)+I_UA(1));
E(1) = x(5);

O(1) = R_D(1);

S(1) = N - E(1) - I_D(1) - I_U(1) - R_D(1) - R_U(1);

F(1) = 0;
C(1) = 0;

%% Calculate time path

for t = 1:T-1;
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
    if Tstart_lock==0;
        lambdat(t) = 0;
        else
            if t < Tstart_lock
            lambdat(t) = 0;
            elseif (t >= Tstart_lock) & (t < Tend_lock);
            lambdat(t) = lambda_P*(1-exp(-eta_L*(t-Tstart_lock)));
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
        betat(t+1) = beta0*a(t)^x(8);
        
end;

a_m = data.a_m(1:length(a));

%dF = F(2:T)-F(1:T-1);
%dF_m = data.dF_m;
O_m = data.O_m;
H_m = data.H_m;
F_m = data.F_m;
%C_D_m = data.C_D_m;


%output = sum( (O_m - O).^2 ) + sum( (H_m - H).^2 ) + sum( (F_m - F).^2 )...
%    + sum( (a_m - a).^2 ) + sum( (dF_m - dF).^2 ) + sum( (C_D_m - C_D).^2 );

output = sum( (O_m - O).^2 ) + 1000*sum( (F_m - F).^2 )...
     +  1000*sum( (a_m - a).^2 );




