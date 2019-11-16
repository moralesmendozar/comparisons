%% UPENN, @Wharton
% Finance 937. 
% Prof. Joao Gomes
% Student: Rodrigo A. Morales M. && Mr. Paw Bednarek
% 
% Based on Jesus Fernandez-Villaverde RBC comparison code
% Okt 5, 2019

% Problem Set 01. Exercise d)

%% 0. Housekeeping

clear all
close all
clc


%%  1. Calibration
clear
clc

rho = 0.8;   % ar(1) parameter of log(a) (productivity)
r = 0.02;    % 1/(1+r) discount rate of firm
bbeta = 1/(1+r);   % Discount factor of the firm
delta = 0.1; % depreciation
theta1 = 0.3; % kapital elasticity (cobbdouglas)
aalpha =theta1;     % Elasticity of output w.r.t. capital
theta2 = 0.6; % labor elasticity (cobbdouglas)
W = 2;   % wage
sigma = 0.1;  % std dev of eps (ar(1)) of log(a) (productivity)
abar = 8; % log(abar) is the mean of log(a), which is ar(1)
b0 = 0;    %parameter0 for cots of investing
b1 = 0.5;   %parameter1 for cots of investing
options = optimset('Display', 'off');  %when solves for kss

% Productivity values
na = 9;  %number of points for logaGrid = [loga0, loga1... loga9]
logabar = log(abar);
deltaLoga = sigma/sqrt(1-rho^2);  %change between each log(a);
%%%%%%%%%%%%%%%%
% For d), only take the meanshock 
vProductivity = logabar;
nGridProductivity = length(vProductivity);

% Transition matrix
if nGridProductivity > 1
    P   = eye(na);
    a_1 = vProductivity(1);
    a_n = vProductivity(na);
    for j = 1:na
        aj = vProductivity(j);
        upperBoundA = (a_1 - rho*aj - (1-rho)*logabar +deltaLoga/2)/sigma;
        P(1,j) = normcdf(upperBoundA);
        lowerboundA = (a_n - rho*aj - (1-rho)*logabar -deltaLoga/2)/sigma;
        P(na,j) = 1-normcdf(lowerboundA);
    end
    for i = 2:(na-1)
        for j = 1:(na)
            ai = vProductivity(i);
            aj = vProductivity(j);
            upperBoundA = (ai - rho*aj - (1-rho)*logabar +deltaLoga/2)/sigma;
            lowerboundA = (ai - rho*aj - (1-rho)*logabar -deltaLoga/2)/sigma;
            P(i,j) = normcdf(upperBoundA)-normcdf(lowerboundA);
        end
    end
    P'
    mTransition   = P;
else
    mTransition = 1;
end


%% 2. Steady State
%define the functions of labor, profit and investment...
labor = @(a,k) (theta1*(k.^theta1)'*a/W).^(1/(1-theta2));
profit = @(a,k) ((k'.^theta1)*a).*(labor(a,k).^theta2) - W*labor(a,k);
investment = @(a,k,kprime) kprime - (1-delta)* k'*ones(size(a,1));
%for d) 1 and 2
phi = @(a,k,kprime) b0 * k'*ones(1,size(a,1)) +...
    b1*(( investment(a,k,kprime)./ ( k'*ones(1,size(a,1)) )- delta).^2).*k'*ones(1,size(a,1));
%for d) (3) ,  b1 =10  if investment is higher than delta*k
phi2 = @(a,k,kprime) b0 * k'*ones(1,size(a,1)) +...
     10*(( investment(a,k,kprime)./ ( k'*ones(1,size(a,1)) )-delta).^2).*k'*ones(1,size(a,1));
%for d) (4)
%phi3 = @(a,k,kprime) b0 * k'*ones(1,size(a,1)) +...
%     b1*(( investment(a,k,kprime)./ ( k'*ones(1,size(a,1)) )-
%     delta).^2).*k'*ones(1,size(a,1));

%find the values of the steady state:
fss =@(a,kss) delta + b0 - (theta1*a*kss^(theta1-1)*labor(a,kss)^theta2)/(1+r);
[capitalSteadyState,fval] = fsolve(@(k)fss(abar,k),30,options);
outputSteadyState = abar*capitalSteadyState^theta1*labor(abar,capitalSteadyState)^theta2;

fprintf(' Output = %2.6f, Capital = %2.6f\n', outputSteadyState, capitalSteadyState); 
fprintf('\n')

% with kss, generate the grid of capital
capMiddle = (aalpha*abar/(r+delta))^(1/(1-aalpha));
kstep = 0.1; %0.00001
%vGridCapital = 0.5*capitalSteadyState:kstep:1.5*capitalSteadyState;
vGridCapital = 0.5*capMiddle:kstep:1.5*capMiddle;
kapitalMax = (abar/(delta))^(1/(1-aalpha));
vGridCapital = (0):kstep:(0.6*kapitalMax);
vGridCapital = (0):kstep:(capitalSteadyState);
vGridCapital = linspace(0.01*capitalSteadyState,capitalSteadyState,200);



nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);

%% Plot of the profit function

figure(1)
plot(vGridCapital, profit(exp(vProductivity(1)),vGridCapital),'-b')
hold on
plot(vGridCapital, profit(exp(vProductivity(nGridProductivity)),vGridCapital),'-r')
hold off
title('profit(k) for lowest & Highest productivity')

%% 3. Required matrices and vectors

mOutput           = zeros(nGridCapital,nGridProductivity);
mValueFunction    = zeros(nGridCapital,nGridProductivity);
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

%% 4. We pre-build output for each point in the grid

profitMatrix = profit(exp(vProductivity),vGridCapital);

%% 5. Main iteration

maxDifference = 10.0;
tolerance = 0.00001;%0.0000001;
iteration = 0;
tic
while (maxDifference>tolerance)  
    
    expectedValueFunction = mValueFunction*mTransition';
    
    for nProductivity = 1:nGridProductivity
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for nCapital = 1:nGridCapital
                        
            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
                
                %consumption = mOutput(nCapital,nProductivity)-vGridCapital(nCapitalNextPeriod);
                prodctvt = exp(vProductivity(nProductivity));
                currentCapital = vGridCapital(nCapital);
                kprimeTomorow = vGridCapital(nCapitalNextPeriod);
                invstm = investment(prodctvt,currentCapital,kprimeTomorow);
                %for 01.(d)(3)
                if invstm > delta*currentCapital
                    dividend = profitMatrix(nCapital,nProductivity) - invstm - phi(prodctvt,currentCapital,kprimeTomorow);
                else
                    dividend = profitMatrix(nCapital,nProductivity) - invstm - phi2(prodctvt,currentCapital,kprimeTomorow);
                end
                valueProvisional = (1-bbeta)*dividend+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);
            
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                    break; % We break when we have achieved the max
                end    
                  
            end
            
            mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar;
            mPolicyFunction(nCapital,nProductivity) = capitalChoice;
            
        end
        
    end
    
    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
    end
           
end

fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
fprintf('\n')
%fprintf(' My check = %2.6f\n', mPolicyFunction(1000,3));
fprintf('\n')

toc

%% 6. Plotting results

figure(2)

subplot(2,1,1)
plot(vGridCapital(1:100),mValueFunction(1:100,:))
xlim([vGridCapital(1) vGridCapital(nGridCapital/2)])
xlabel('k')
title('Value Function')

subplot(2,1,2)
plot(vGridCapital(1:100),mPolicyFunction(1:100,:))
xlim([vGridCapital(1) vGridCapital(nGridCapital/2)])
xlabel('k')
title('Policy Function')

% vExactPolicyFunction = aalpha*bbeta.*(vGridCapital.^aalpha);
% 
% subplot(3,1,3)
% plot((100.*(vExactPolicyFunction'-mPolicyFunction(:,1))./mPolicyFunction(:,1)))
% title('Comparison of Exact and Approximated Policy Function')

%set(gcf,'PaperOrientation','landscape','PaperPosition',...
%[-0.9 -0.5 12.75 9])
%print('-dpdf','Figure1.pdf')

%% Plotting optimal investment and financing for lowest, avg and highest
mOptimalInvestment = zeros(nGridProductivity, nGridCapital);
mFinancing = zeros(nGridProductivity, nGridCapital);
for nProductivity = 1:nGridProductivity
    for nCapital = 1:nGridCapital
        prodctvt = exp(vProductivity(nProductivity));
        currentCapital = vGridCapital(nCapital);
        captomorrow = mPolicyFunction(nCapital,nProductivity);
        invstmnt = investment(prodctvt,currentCapital,captomorrow);
        mOptimalInvestment(nCapital,nProductivity) = invstmnt;
        profits = profitMatrix(nCapital,nProductivity);
        costOfInvstmnt = phi(prodctvt,currentCapital,captomorrow);
        dividend = profits - invstmnt - costOfInvstmnt;
        mFinancing(nCapital,nProductivity) = max(-dividend,0);
    end
end
%shock a
figure(3)
% subplot(2,1,1)
plot(vGridCapital(1:100),mOptimalInvestment(1:100,1))
xlabel('k')
title('d)01) Optimal Invesment for the mean shock only')

% hold on
% plot(vGridCapital,mOptimalInvestment(:,5))
% plot(vGridCapital,mOptimalInvestment(:,1))
% legend('high a', 'middle a', 'low a')
% hold off

% subplot(2,1,2)
% plot(vGridCapital(1:100),mFinancing(1:100,1))
% xlabel('k')
% title('d) Financing (-d(a,k) when positive) for the mean shock only')

% hold on
% plot(vGridCapital,mFinancing(:,5))
% plot(vGridCapital,mFinancing(:,1))
% legend('high a', 'middle a', 'low a')
% hold off