clear all; clf; %We clear the workspace

tEnd = 60*120; %tEnd is the time at which the simulation stops
tStep = 0.001;
j=1;


for i=[0.1 1 3 5 7]
    
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

% values on wants to change without changing the Alpha' and without
% affecting Beta'
alphaM = 0.5*i;    
Km = 40*i;

gammaM = 1/120;
gammaP = 1/600;
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

if i==0.1
    Species01=Species; 
    t01=t;
elseif i==1
   Species1=Species; 
   t1=t;
elseif i==3
   Species3=Species; 
   t3=t;
elseif i==5
   Species5=Species; 
   t5=t;
elseif i==7
   Species7=Species;
   t7=t;
else 
    disp('ERROR')
end


disp(j)
j=j+1;

end

%% plotting

subplot(5,1,1)
hold on
plot(t01,Species01(:,4))
plot(t01,Species01(:,5))
plot(t01,Species01(:,6))
xlabel('Time')
ylabel('Concentrations')
legend('LacI Protein','TetR Protein','CI Protein')
title('for i = 0.1')


subplot(5,1,2)
hold on
plot(t1,Species1(:,4))
plot(t1,Species1(:,5))
plot(t1,Species1(:,6))
xlabel('Time')
ylabel('Concentrations')
legend('LacI Protein','TetR Protein','CI Protein')
title('for i = 1')

subplot(5,1,3)
hold on
plot(t3,Species3(:,4))
plot(t3,Species3(:,5))
plot(t3,Species3(:,6))
xlabel('Time')
ylabel('Concentrations')
legend('LacI Protein','TetR Protein','CI Protein')
title('for i = 3')

subplot(5,1,4)
hold on
plot(t5,Species5(:,4))
plot(t5,Species5(:,5))
plot(t5,Species5(:,6))
xlabel('Time')
ylabel('Concentrations')
legend('LacI Protein','TetR Protein','CI Protein')
title('for i = 5')

subplot(5,1,5)
hold on
plot(t7,Species7(:,4))
plot(t7,Species7(:,5))
plot(t7,Species7(:,6))
xlabel('Time')
ylabel('Concentrations')
legend('LacI Protein','TetR Protein','CI Protein')
title('for i = 7')



