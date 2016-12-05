% REPRESSILATOR
clc;
clear all;

%% Parameters

alpha = [10^0 10^(0.5) 10^1 10^(1.5) 10^2 10^(2.5) 10^3 10^(3.5) 10^4 10^(4.5)];% Promoter strength, from 5,10^4 to 0,5
alphav = [0 0.5 1 1.5 2 2.5 3 3.5 4 4.5]; %to use in the plot
beta = [10^0 10^(0.5) 10^1 10^(1.5) 10^2 10^(2.5) 10^3];% Average translation efficiency (IN PROT/TRANSCRIPT)
betav = [0 0.5 1 1.5 2 2.5 3]; %to use in the plot
alpha0 = 0 ;% Leakiness of the promoter. can be 0, 10^(-3), ...
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

Oscil=zeros(2,1);
Noscil=zeros(2,1);

%% Modifying alpha and beta

for j=1:length(alpha)
    
    alphaj=alpha(j);
    
    for k=1:length(beta)
        
        betaj=beta(k);
        
        %% Reseting the experiment
        t=0;
        tEnd = 10;
        tStep = 0.001;
        
        p(1) = PlacI;
        p(2) = PtetR;
        p(3) = PcI;
        m(1) = MlacI;
        m(2) = MtetR;
        m(3) = McI;
        
        P = zeros(tEnd/tStep, 3); % Number of proteins
        M = zeros(tEnd/tStep, 3); % Number of mRNA
        
        %% Reactions
        
        i=1;
        
        while t <= tEnd
            P(i, 1) = p(1);
            P(i, 2) = p(2);
            P(i, 3) = p(3);
            M(i, 1) = m(1);
            M(i, 2) = m(2);
            M(i, 3) = m(3); % storage of the values along the experiment
            
            
            %%% mRNA
            m(1) = (- m(1) + alphaj./(1 + p(2).^n) + alpha0)*tStep + m(1);
            m(2) = (- m(2) + alphaj./(1 + p(3).^n) + alpha0)*tStep + m(2);
            m(3) = (- m(3) + alphaj./(1 + p(1).^n) + alpha0)*tStep + m(3);
            
            %%% Proteins
            p(1) = (- betaj*(p(1)-m(1)))*tStep + p(1);
            p(2) = (- betaj*(p(2)-m(2)))*tStep + p(2);
            p(3) = (- betaj*(p(3)-m(3)))*tStep + p(3);
            
            i = i + 1;
            t(end +1) = t(end) + tStep;
        end
        
        
        M(end+1,:) = m(end,:);
        P(end+1,:) = p(end,:); % last values
        
        x=std(P(8000:end,1))/mean(P(8000:end,1));
        
        OSCindex(j,k)=x;
        
        %% Dividing the values in two categories : the ones with oscillations, the ones without
        
        if OSCindex(j,k)>0.1
            Oscil(:,end+1)=[alphav(j);betav(k)];
                        
        else 
            Noscil(:,end+1)=[alphav(j);betav(k)];
      
        end
        
        hold on
        plot(Oscil(1,:),Oscil(2,:),'ro')
        plot(Noscil(1,:),Noscil(2,:),'b*')
        xlabel('Order of magnitude of Alpha (10^x)')
        ylabel('Order of magnitude of Beta (10^y)')
        legend('Oscillations','No oscillation') 
        title('Phase diagram for parameters Alpha and Beta')
        
        
    end
end

%% Plots
