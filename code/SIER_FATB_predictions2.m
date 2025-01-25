




function output = SIER_FATB_predictions2(para);



%% Step 4: display estimation results (ans save more results)

sigma = para.sigma;  %rate at which exposed become symptomatic
gamma = para.gamma;    %rate at which infectious become "resolving"
theta = para.theta;   %rate at which disease is resolved (transit to recovered or dead)

delta0 = para.delta0;   %death rate 0.08=0.8% death rate
delta1 = para.delta1;
eta_d = para.eta_d;

psi = para.psi;
R_0 = para.R_0;        %initial virality
beta0=R_0*gamma;         % initial beta
alpha = para.alpha;    %Probability of NOT displaying symptoms

N=para.N; %Population
T=para.T; %Time horizon

Tfor = para.Tfor;
Tpolicy=para.Tpolicy;
%% Policy parameters

    %Policy 1: Testing & quarantining
        tauR = para.tauR;      %random testing rate
        tauS1 = para.tauS1;    %testing of symptomatic infectious people
        tauS2 = para.tauS2;    %testing of symptomatic infectious people        

        %Tstart_test = 60;
        %kappat = kappa.*[zeros(Tstart_test,1); ones(T-Tstart_test,1)];
        
    %Policy 2: Forced lockdown
        lambda_P=para.lambda_P; %forced lockdown, e.g. lambda=0.5 locks down 40% of the workforce    
        Tstart_lock = para.Tstart_lock;
        Tend_lock = para.Tend_lock;
        eta_L = para.eta_L;
       
   %Policy 3: Information policies
        a0 = para.a0;      %constant
        a1 = para.a1;     %confirmed cases
        a2 = para.a2;    %hospitalized
        a3 = para.a3;      %deaths
        

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
deltat = zeros(T,1);

ls = S;
betat = S;
    betat(1) = beta0;
a(1) = a0;

%% Initialization of disease
R_DS(1) = 1;                          %one initial case (hospitalized)
R_DA(1) = 0;                          %one initial case
    R_D(1) = R_DS(1)+R_DA(1);

R_US(1) = para.R_US_1;    %implied non-hospitalized resolving cases
R_UA(1) = para.R_UA_1;
    R_U(1) = R_US(1)+R_UA(1);

    R_S(1) =R_DS(1)+R_US(1);
    R_A(1) = R_UA(1)+R_DA(1);
    
I_DS(1) = 0;
I_DA(1) = 0;
    I_D(1) = I_DS(1)+I_DA(1);

I_US(1) = para.I_US_1;
I_UA(1) = para.I_UA_1;
    I_U(1) = I_US(1)+I_UA(1);

    I_S(1) = I_DS(1)+I_US(1);
    I_A(1) = I_DA(1)+I_UA(1);
    
E_D(1) = 0;
%E(1) = (1/sigma)*(I_DS(1)+I_DA(1)+I_US(1)+I_UA(1));
E(1) = para.E_1;

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
        deltat(t) =delta0*exp(-eta_d*(t)) +delta1*(1-exp(-eta_d*(t)));
        F(t+1) = F(t) + deltat(t)*(theta*R_US(t)+theta*R_DS(t));
        
        %ReCovered
        C(t+1) = C(t) + theta*R_UA(t) + (1-deltat(t))*theta*R_US(t);
        C_D(t+1) = C_D(t) + theta*R_DA(t) + (1-deltat(t))*theta*R_DS(t);

        
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
            lambdat(t) = para.lambda_P*(1-exp(-eta_L*(t-Tstart_lock)));
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
 


        % Avoidance health effect
        betat(t+1) = beta0*a(t)^psi;
        
end 


output(:,1) = O;
output(:,2) = F;
output(:,3) = a;

end


