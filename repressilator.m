% REPRESSILATOR

clear all; 

%% Parameters

alpha = 0.5 ;% Promoter strength, from 5,10^4 to 0,5
alpha0 = 10^2 ;% Leakiness of the promoter. can be 0, 10^(-3), ...
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
MlacI = 1;
MtetR = 2;
McI = 3;

% Denominations 
p(1) = PlacI;
p(2) = PtetR;
p(3) = PcI;
m(1) = MlacI;
m(2) = MtetR;
m(3) = McI;

t=0;
tEnd = 10;
tStep = 0.1;

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
    dmLac = - m(1) + alpha./(1 + p(2).^n) + alpha0;
    dmTet = - m(2) + alpha./(1 + p(3).^n) + alpha0;
    dmI = - m(3) + alpha./(1 + p(1).^n) + alpha0;
    
    %%% Proteins
    dpLac = - beta*(p(1)-m(1));
    dpTet = - beta*(p(2)-m(2));
    dpI = - beta*(p(3)-m(3));
    
    
    m(1) = dmLac;
    m(2) = dmTet;
    m(3) = dmI;
    p(1) = dpLac;
    p(2) = dpTet;
    p(3) = dpI;
    
    i = i + 1;
    t(end +1) = t(end) + tStep; 
end
M(end+1,:) = m(end,:);
P(end+1,:) = p(end,:);


subplot(2,1,1)
plot(t,M(:,1), 'c')
hold on;
plot(t,M(:,2), 'g')
hold on;
plot(t,M(:,3))
%
%subplot(2,1,2)
%plot(t,P(:,1))
%hold on;
%plot(t,P(:,2))
%plot(t,P(:,3))

pause;