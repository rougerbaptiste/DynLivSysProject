clear all; clf; %We clear the workspace


tEnd = 60*200; %tEnd is the time at which the simulation stops
tStep = 0.001;

% Here we have to put the initial concentrations of the species
% Proteins concentrations
PlacI = 100;
PtetR = 100;
PcI = 100;

stoich =  [ 1 0 0 0 0 0;
            0 1 0 0 0 0;
            0 0 1 0 0 0;
            -1 0 0 0 0 0;
            0 -1 0 0 0 0;
            0 0 -1 0 0 0;
            0 0 0 1 0 0;
            0 0 0 0 1 0;
            0 0 0 0 0 1;
            0 0 0 -1 0 0;
            0 0 0 0 -1 0;
            0 0 0 0 0 -1;
];

% mRNAs concentrations
MlacI = 0;
MtetR = 0;
McI = 0;

% Denominations
species(1) = MlacI;
species(2) = MtetR;
species(3) = McI;
species(4) = PlacI;
species(5) = PtetR;
species(6) = PcI;



gammaM = 1/120;
gammaP = 1/600;
Km = 40;
alphaM = 0.5;
alphaP = 1/6;
leakiness = 5e-04;
n = 2; %hill coefficient

Species = species; %will contain the concentrations of the species over the time


t = 0; % t is the vector containing the times
while t(end) <= tEnd

    propMlacIProd = (alphaM)/(1 + (species(5) / Km)^n ) + leakiness;
    propMtetRProd = (alphaM)/(1 + (species(6) / Km)^n ) + leakiness;
    propMcIProd = (alphaM)/(1 + (species(4) / Km)^n ) + leakiness;
    propMlacIDest = gammaM*species(1);
    propMtetRDest = gammaM*species(2);
    propMcIDest = gammaM*species(3);

    propPlacIProd = alphaP*species(1);
    propPtetRProd = alphaP*species(2);
    propPcIProd = alphaP*species(3);
    propPlacIDest = gammaM*species(4);
    propPtetRDest = gammaM*species(5);
    propPcIDest = gammaM*species(6);

    props = [propMlacIProd propMtetRProd propMcIProd propMlacIDest propMtetRDest propMcIDest propPlacIProd propPtetRProd propPcIProd propPlacIDest propPtetRDest propPcIDest];

    propsum = cumsum(props);
    pickingVector = propsum/propsum(end);

    rand1 = rand();%rand for reaction

    selectedReact = find(pickingVector>rand1,1);

    species = species + stoich(selectedReact,:);


    Species(end+1,:) = species;





    rand2 = rand();% rand for time




    deltaT = - log(rand2)./sum(props);
    t(end+1) = t(end) + deltaT;

end


%% Autocorrelation

Delta=zeros();
MeanLacI=zeros();
C=zeros();
tend=10;
tstep=0.5; %for autocorrelation
Autocorr=zeros(1,tend/tstep);
j=1;

MeanLacI=mean(Species(:,4)); 

for delta = 300*tstep:tstep:length(Species)  
    
    if 1+j<=length(Species)
        
        C= mean((Species(1+j,4)-MeanLacI) * (Species(1,4)-MeanLacI) )/ (MeanLacI)^2;
        
        Autocorr(j)=C;
        Delta(end+1)=delta;
    end
    
    j=j+1;
end


plot(Delta(2:end-1),Autocorr(2:end))