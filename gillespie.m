clear all; clf; %We clear the workspace

tEnd = 10; %tEnd is the time at which the simulation stops
tStep = 0.001

% Here we have to put the initial concentrations of the species


t = 0; % t is the vector containing the times
while t(end) <= tEnd

    rand1 = rand();%rand for reaction
    rand2 = rand();% rand for time

    deltaT = -log(rand2)
    t(end+1) = t(end) + deltaT
end
