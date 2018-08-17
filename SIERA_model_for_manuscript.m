%% repeated dosing of endectocide and bednet delivery
% * SEIRA_humans_livestock_periods_search.m introduces two control mechanisms for mosquitos.
% One delivers N doses of a drug to humans and/or cattle at pre-determined
% intervals, the other is the use of LLINs. Output is time courses of
% malaria prevalence in humans and mosquito populations.
% * Humans are represented by 5 categories: susceptible, exposed, infectious,
% immune, and asymptomatic.
% * Mosquitoes are represented by 3 categories: susceptible, exposed, and
% infectious
% * The time search is finer than previous iterations where levels of week,
% month, 1/4 annual, 1/2 annual and annual were implemented. This is to try
% to produce smoother contour maps for the dosing frequency needed to
% achieve target relative reduction in malaria prevalence. 
% * Biting is no longer a function of alpha and beta. Instead it is linearly
% correlated to human host prevalence.

% by Hannah Meredith
% last updated: August 105, 2018

tic
%% Parameters that change based on treatment

N_death =0.73; % LD50 of standard susceptible mosquitoes to LLIN [%]
hill_N = 2; % hill coefficient for pyrethroids: 3 for susceptible, 2 for resistant mosquitoes
hill_D = 2; % hill coefficient for systemic insecticides
hl_n = 1906; % half life net (days)
mu_d = 0.6; % mosquito death by insecticide [0.6 day ^-1]
mu_n = 0.1; % mosquito death rate by bed net [0.1 - 0.5 day^1]

% Parameters for control methods
t_ss = 1000;        % time malaria prevalence reaches steady state
period_N = 1095;    % recommended time between net delivery is 3 years (days)
net = 2;         % concentration of permethrin in new LLINs [%]

% initial density for history function
E_h0 = 0.01; % exposed human
I_h0 = 0.01; % infected human
R_h0 = 0.01; % recovered human
A_h0 = 0.01; % asymptomatic human
E_m0 = 0.01; % exposed mosquito
I_m0 = 0.01; % infected mosquito
D_l0 = 0; % drug concentration in livestock
D_h0 = 0; % drug concentration in human
N0 = 0; % drug concentraiotn in net
y0 =[E_h0; I_h0; R_h0; A_h0; E_m0; I_m0; D_l0; D_h0; N0];

%% simulation
cmax_range = logspace(0,4,3); %logspace(0,4,10); %linspace(1,3000,10) % [ng/mL]
period_D_range = [365/52, linspace(365/12,365,3)];%[365/52, linspace(365/12,365,12)];%[365/52,365/12,365/4,365/2,365];%linspace(7, 365, 15); % establish range of dosing periods [days]
halflife_d_range= logspace(-1,2, 3); %logspace(-1,2, 10); %[0.5, 2.5, 7.5, 12.5, 25];%linspace(0.3,40, 5); % range of 1/2 lifes [days]
D_death_range=[10, 25, 100, 500, 1100];%linspace(1,3000,10); % range of LC50s
m_range = [5, 20, 40];
p_H_range = [0.35, 0.8];
toggle1=0; %switches off at preset malprev level
target=10;% target reduction in malaria prevalence rate


C_n = 1 * 0.75; % p21. coverege by endectocides (0:1) * max efficacy of newly applied net (0.75)
C_l = 1 * 1; % p22. coverege of livestock by endectocides (0:1) * max efficacy of newly applied drug (1)
C_h = 0 * 1; % p22. coverege of humans by endectocides (0:1) * max efficacy of newly applied drug (1)

for h = 1:length(D_death_range)%1:length(D_death_range) % threshold concentration of drug for mosquito death
    D_death = D_death_range(h);
    m = 20;%m_range(h);
    p_H = 0.5;%p_H_range(h);

    periodUse = zeros(length(halflife_d_range), length(cmax_range));
    rel_mal_H = zeros(length(halflife_d_range), length(cmax_range), length(period_D_range));
    
    for j = 1:length(halflife_d_range) % half-life of drugs
        hl_d = halflife_d_range(j);%
        
        for k = 1:length(cmax_range) % max concentration of drugs (initial concentration in  plasma)
            dose = cmax_range(k);%
            
            
            for i = length(period_D_range):-1:1   % time between doses
                period_D = period_D_range(i); %% 
                
                
                p = [hl_d; hl_n; mu_d; mu_n; D_death; N_death; C_n; C_l; C_h; p_H; m; hill_N; hill_D];
                control = [t_ss, period_D, period_N, dose, net];
                
                [T_base,I_base, base_ss]=baseline(y0,p,control);
                [T_net,I_net,mal_prev_net]=nets(y0,p,control);
                frax_net = I_net./base_ss;
%                 frax_net = I_net./I_base;
                
                [T,Y,doses] = dosing(y0, p, control);
                
                E_h = Y(:,1);
                I_h = Y(:,2);
                R_h = Y(:,3);
                A_h = Y(:,4);
                E_m = Y(:,5);
                I_m = Y(:,6);
                D_l = Y(:,7);
                D_h = Y(:,8);
                N = Y(:,9);
                time = T;
                T_h = I_h + A_h;
                frax_si = T_h/base_ss;  % fraction of malaria prevalence with nets + systemic insecticide present compared to no control
%                 frax_si = T_h./I_base;  %
                S_h = 1 - (E_h + I_h + R_h + A_h);
                
                figure(1)
                hold on
                plot(time, T_h, 'LineWidth',2)
%                 plot (time, S_h, time, E_h, time, I_h, time, R_h,time, A_h)
%                 subplot(1,2,1)
%                 hold on
%                 plot (T_net, frax_net,time, frax_si,'LineWidth', 3)%  plot (time, frax_si,'LineWidth', 3)
                xlim([1000 4650])
                ylim([0 1])
                xticks([1000+0 1000+5*365 1000+10*365 1000+15*365])
                xticklabels({'0','5','10','15'})
%                 xlabel('Time (years)')
%                 ylabel('Malaria prevalence (%)')
%                 axis square
               
                mal_H(j,k) = T_h(end); %mal_H(j,k,i) = T_h(end);
                rel_mal_H(j,k,i)= (mal_prev_net-T_h(end))./(mal_prev_net)*100;%rel_mal_H(j,k,i)= (mal_prev_net-T_h(end))./(mal_prev_net)*100;
                rel_mal_H(rel_mal_H<0)=0;
                                
                [i,k,j,h]
                                
                if rel_mal_H(j,k,i) >= target
                    periodUse(j,k)=period_D;
                    break
                end
            end
        end        
    end  
    
%     
        figure(1)
        hold on
        subplot(2,3,2)
        v=[365/52,365/12,365/4,365/2,365];%period_D_range;
        contour(halflife_d_range, cmax_range, periodUse,v)
        shading flat
        colormap('parula')
        axis square
        
%         figure(8)
%         hold on
%         contour(halflife_d_range, cmax_range, rel_mal_H)%contour(halflife_d_range, cmax_range, matrix)
%         shading flat
%         colormap('parula')
%         axis square
end

toc

%% functions

% baseline simulation

function [T, I_base, base_ss] = baseline(y0, p, c)
p(7)=0;
p(8)=0;
p(9)=0;

[T,Y] = dosing(y0, p, c);
I_base = Y(:,2)+Y(:,4);
base_ss = I_base(end);
end

% nets alone simulation

function[T, I_net, mal_prev_net] = nets(y0, p, c)
p(7)=0.75;
p(8)=0;
p(9)=0;
[T,Y] = dosing(y0, p, c);
I_net = Y(:,2)+Y(:,4);
mal_prev_net= I_net(end);
end

% Dosing function - solves time course in piece-meal fashion for multi-dosing
function [overallT, overallY, doseNum, netNum] = dosing(y0, p, c)

t_ss = c(1);
period_D = c(2);
period_N = c(3);
dose = c(4);
net = c(5);

hl_d = p(1);
hl_n = p(2);

% Allow time to pass while malaria prevalence reaches steady state; no control methods in place
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t, y] = ode45(@transmission, [0 t_ss], y0, options, p);
overallY = y(:,:);
overallT = t;

%first dose for both campaigns. Assumes drug will be dosed again first

t_startD = overallT(end,1);
t_startN = overallT(end,1);
t_endD = overallT(end,1) + period_D;
t_endN = overallT(end,1) + period_N;
y0 = [y(end,1);y(end,2);y(end,3);y(end,4);y(end,5);y(end,6);y(end,7)+ dose;y(end,8) + dose; net];
[t, y] = ode45(@transmission, [t_startD t_endD], y0, options, p);

overallY = [overallY; y(2:end, :)];
overallT = [overallT; t(2:end, :)];
t_startD = overallT(end);
t_endD = overallT(end) + period_D;
t_startN = t_startN + period_N;
t_endN = t_startN + 2*period_N;
doseNum = 1;
netNum = 1;
treat_pd = t_ss + 10*365;

% Prevalence has reached steady state. Begin dosing regimen
while (overallT(end) < treat_pd) % run simulation for X years or until malaria is below threshold
    
    if (t_endD > treat_pd || t_endD > treat_pd)
        t_endD = treat_pd;
        t_endN = treat_pd;
    end
    
    
    if (abs(t_startD - t_startN)<0.0001)
        
        y0 = [y(end,1);y(end,2);y(end,3);y(end,4);y(end,5);y(end,6);y(end,7)+ dose; y(end,8)+ dose; net];
        [t, y] = ode45(@transmission, [t_startD t_endD], y0, options, p);
        
        overallY = [overallY; y(2:end, :)];
        overallT = [overallT; t(2:end, :)];
        t_startD = overallT(end);
        t_endD = overallT(end) + period_D;
        t_startN = t_startN + period_N;
        t_endN = t_startN + period_N;
        netNum = netNum + 1;
        doseNum = doseNum + 1;
        
        
    else if (t_startD < t_startN && t_endD < t_startN)
            
            y0 = [y(end,1);y(end,2);y(end,3);y(end,4);y(end,5);y(end,6);y(end,7)+dose;y(end,8)+ dose;y(end,9)];
            [t, y] = ode45(@transmission, [t_startD t_endD], y0, options, p);
            
            overallY = [overallY; y(2:end, :)];
            overallT = [overallT; t(2:end, :)];
            doseNum = doseNum + 1;
            t_startD = overallT(end);
            t_endD = overallT(end) + period_D;
            
        else if (t_startD < t_startN && t_endD >= t_startN)
                
                y0 = [y(end,1);y(end,2);y(end,3);y(end,4);y(end,5);y(end,6);y(end,7)+dose;y(end,8)+dose;y(end,9)];
                [t, y] = ode45(@transmission, [t_startD t_startN], y0, options, p);
                
                overallY = [overallY; y(2:end, :)];
                overallT = [overallT; t(2:end, :)];
                doseNum = doseNum + 1;
                t_startD = t_endD;
                t_endD = t_endD + period_D;
                t_startN = overallT(end);
                t_endN = overallT(end) + period_N;
                
            else if (t_startN < t_startD)
                    y0 = [y(end,1); y(end,2); y(end,3); y(end,4); y(end,5); y(end,6); y(end,7); y(end,8); net];
                    [t, y] = ode45(@transmission, [t_startN t_startD], y0, options, p);
                    
                    overallY = [overallY; y(2:end, :)];
                    overallT = [overallT; t(2:end, :)];
                    netNum = netNum + 1;
                    t_startN = t_startN + period_N;
                    t_endN = t_startN + period_N;
                    
                    y0 = [y(end,1); y(end,2); y(end,3); y(end,4); y(end,5); y(end,6); y(end,7)+dose; y(end,8)+dose; y(end,9)];
                    [t, y] = ode45(@transmission, [t_startD t_endD], y0, options, p);
                    
                    overallY = [overallY; y(2:end, :)];
                    overallT = [overallT; t(2:end, :)];
                    doseNum = doseNum +1;
                    t_startD = overallT(end);
                    t_endD = overallT(end) + period_D;                    
                    
                end
            end
        end
    end
end

% figure()
% hold on
% Y=overallY(:,2)+overallY(:,4);
% subplot(1,2,1)
% plot(overallT, Y)
% subplot(1,2,2)
% hold on
% yyaxis left
% plot(overallT,overallY(:,7),'LineWidth', 1)
% yyaxis right
% plot(overallT, overallY(:,9), 'LineWidth', 1)

end

% function to solve malaria transmission
function dydt = transmission(t, Y, p)
E_h = Y(1);
I_h = Y(2);
R_h = Y(3);
A_h = Y(4);
E_m = Y(5);
I_m = Y(6);
D_l = Y(7);
D_h = Y(8);
N = Y(9);

% define parameters
a = 0.2; % biting frequency [ 0.01 – 0.5 day-1]
b = 0.5; % proportion of bites that produce infection in humans [0.2 – 0.5]
c = 0.5; % proportion of bites that produce infection in mosquitoes [0.5]
r = 0.01; % rate of recovery [0.005 – 0.05 day^-1]
mu_1= 1/21900; % net change in human population [1/21900 day^-1]
mu_2 = 0.12; % birth (maturation)/death rate of mosquitoes [0.05 - 0.5 day^-1]
tau_m = 10; % latent period for mosquitoes [5 - 15 days]
tau_h = 21; % latend period for humans [10 - 100 days]
q1 = 1/200; % rate of gaining immunity
q2 = 1/1000; % rate of losing immunity
sigma = 0.25; % adjustment factor for asymp transmissibility to vector
theta = 0.5; % level of reduced susceptibility to secondary infection

hl_d = p(1); % half life of drug
hl_n = p(2); % hald life of net
k_d = log(2)/hl_d; % elimination rate constant of drug (days ^-1)
k_n = log(2)/hl_n; % decay rate constant of LLIN [days^-1]
mu_d = p(3); % mosquito death by insecticide [0.6 day ^-1]
mu_n = p(4); % mosquito death by LLIN 
D_death = p(5);% LC50 of drug
N_death = p(6); % LC50 of LLIN
C1 = p(7);% coverege by nets (0:1) * max efficacy of newly applied drug (0.75)
C2 = p(8);% livestock coverege by endectocides (0:1) * max efficacy of newly applied drug (1)
C3 = p(9);% human coverege by endectocides (0:1) * max efficacy of newly applied drug (1)
p_h = p(10); % proportion of bites on humans
m = p(11); % ratio of mosquitoes to blood meal sources [0.5 - 40]
hill_N = p(12);
hill_D = p(13);

b_h = a * b * (1 - C1 * N.^hill_N ./(N.^hill_N + N_death^hill_N)); % diminised tramsmission potential from vector to human due to LLIN reducing bite rate
b_m = a * c * (1 - C1 * N.^hill_N ./(N.^hill_N + N_death^hill_N)); % diminised tramsmission potential from human to vector due to LLIN reducing bite rate
mu_2c = mu_2 + 1/3 * (p_h * C1 * N.^hill_N /(N.^hill_N + N_death.^hill_N)*mu_n +(1 - p_h) * C2 * D_l.^hill_D/(D_l.^hill_D + D_death.^hill_D)*mu_d+(p_h) * C3 * D_h.^hill_D/(D_h.^hill_D + D_death.^hill_D)*mu_d); % effect of drug on mosquito death rate

S_h = 1 - E_h - I_h - R_h - A_h;
S_m = 1 - E_m - I_m;

dE_h = m * b_h * p_h * S_h * I_m  - (1/tau_h + r + mu_1)* E_h;
dI_h = 1/tau_h * E_h - (r + q1 + mu_1) * I_h;
dR_h = q1 * (I_h + A_h) - (theta * b_h * m * p_h * I_m + q2 + mu_1) * R_h;
dA_h = theta * b_h * m * p_h * I_m * R_h - (q2 + r + mu_1)*A_h; % change back to q1
 
dE_m = b_m * p_h * (I_h +sigma * A_h) * S_m - (1/tau_m + mu_2c) * E_m;
dI_m = (1/tau_m) * E_m - mu_2c * I_m;
dD_l = - k_d*D_l;
dD_h = - k_d*D_h;
dN = - k_n*N;

dydt = [dE_h; dI_h; dR_h; dA_h; dE_m; dI_m; dD_l;dD_h; dN];
end
