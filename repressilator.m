% REPRESSILATOR
clf;
clear all; 
ducon
%% Parameters

alpha = 10^2 ;% Promoter strength, from 5,10^4 to 0,5
alpha0 = 0 ;% Leakiness of the promoter. can be 0, 10^(-3), ...
beta = 10^2 ;% Average translation efficiency (IN PROT/TRANSCRIPT)
n = 2; % Hill coefficient
% Km = 40 ; % (IN MONOMERS/CELL)
% Phl = 10 % protein half life (IN MIN)
% Mhl = 2 %mRNA half life (IN MIN)

%% Variables

% Proteins concentrations
PlacI = 0;
PtetR = 0;
PcI = 0;

% mRNAs concentrations 
MlacI = 10;
MtetR = 20;
McI = 30;

% Denominations 
p(1) = PlacI;
p(2) = PtetR;
p(3) = PcI;
m(1) = MlacI;
m(2) = MtetR;
m(3) = McI;

t=0;
tEnd = 10;
tStep = 0.001;

P = zeros(tEnd/tStep, 3);
M = zeros(tEnd/tStep, 3);

i=1;
while t <= tEnd
    P(i, 1) = p(1);
    P(i, 2) = p(2);
    P(i, 3) = p(3);
    M(i, 1) = m(1);
    M(i, 2) = m(2);
    M(i, 3) = m(3);
        
    
    %%% mRNA
    m(1) = (- m(1) + alpha./(1 + p(2).^n) + alpha0)*tStep + m(1);
    m(2) = (- m(2) + alpha./(1 + p(3).^n) + alpha0)*tStep + m(2);
    m(3) = (- m(3) + alpha./(1 + p(1).^n) + alpha0)*tStep + m(3);
    
    %%% Proteins
    p(1) = (- beta*(p(1)-m(1)))*tStep + p(1);
    p(2) = (- beta*(p(2)-m(2)))*tStep + p(2);
    p(3) = (- beta*(p(3)-m(3)))*tStep + p(3);
    
    i = i + 1;
    t(end +1) = t(end) + tStep; 
end
M(end+1,:) = m(end,:);
P(end+1,:) = p(end,:);


subplot(2,1,1)
plot(t,M(:,1))
hold on;
plot(t,M(:,2))
plot(t,M(:,3))

subplot(2,1,2)
plot(t,P(:,1))
hold on;
plot(t,P(:,2))
plot(t,P(:,3))