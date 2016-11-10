% REPRESSILATOR

clearall; 

%% Parameters

alpha = ;% Promoter strength, from 5,10^4 to 0,5
alphai = 0 ;% Leakiness of the promoter. can be 0, 10^(-3), ...
beta = ;% Average translation efficiency (IN PROT/TRANSCRIPT)
n = 2; % Hill coefficient
% Km = 40 ; % (IN MONOMERS/CELL)
% Phl = 10 % protein half life (IN MIN)
% Mhl = 2 %mRNA half life (IN MIN)

%% Variables

% Proteins concentrations
PlacI = ;
PtetR = ;
PcI = ;

% mRNAs concentrations 
MlacI = ;
MtetR = ;
McI = ;

% Denominations 
p(1) = PlacI;
p(2) = PtetR;
p(3) = PcI;
m(1) = MlacI;
m(2) = MtetR;
m(3) = McI; 

