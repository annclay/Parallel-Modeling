function [t, population] = parallell_chikungunya_dengue_model(MaxTime,initialStates,parallelOn, r,TranChik,TranDen,recChik,recDen,nu,mu, pSeqChik, pSeqRecovery, symptChik, symptDen0, symptDen1, HRChik, HRDen,limRepeatInfDen)

%Chikungunya Dengue Parallel Transmission models
%Anneke Claypool
%Version 2
%June 27, 2020

% This code runs the differential equations for the parallel chikungunya
% and dengue models, incorporating competing mortality. 

% Vaccination is not considered in this code. 

% This model is adapted from the MATLAB version 4.2 from page 123 of 
% "Modeling Infectious Disease in humans and animals" 
% by Keeling & Rohani.
% 
%% Parallel Modeling- Competing mortality on
compMortOn = parallelOn;

%% Combined Model
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nsteps= MaxTime;
tspan = linspace(1,MaxTime, nsteps); %create timespan with results for each time step
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

[t, population]=ode23t(@Diff,tspan,initialStates,options); % Use ode23t, The main ODE, Originally ode23tb

%% ODE Functions

% Calculates the differential rates used in the integration for chikungunya and dengue.
    function dPop=Diff(t, pop)    

        %% Chikungunya model
        %dS_Hchik/dt 
        S_Hchik = pop(1); 
        %I_H0chik
        I_H0chik = pop(2); 
        %dI_H1chik/dt
        I_H1chik = pop(3);  
        %dR_H0chik/dt
        R_H0chik = pop(4); 
        %dR_H1chik/dt
        R_H1chik = pop(5);
         
        %% Dengue Model
        %S_H0den/dt 
        S_H0den = pop(6);
        %S_H1den/dt
        S_H1den = pop(7);
        %I_H0den/dt
        I_H0den = pop(8);
        %I_H1den/dt
        I_H1den = pop(9);

        %% Chikungunya Transmission Mosquitos
        %dS_Mchik/dt
        S_Mchik = pop(10);
        %dI_Mchik/dt
        I_Mchik = pop(11);
        
        %% Dengue Transmission Mosquitos
        %dS_Mden/dt
        S_Mden = pop(12);
        %I_Mden/dt
        I_Mden = pop(13);
        
        % %Total Mosquito Population
         TotalMosPopChik = sum(pop(10:11));
         TotalMosPopDen = sum(pop(12:13));

        %% Extra data
        %deathsTotChikModel
         deathsTotChik = pop(14);
        %deathsTotDenModel
         deathsTotDen = pop(15); 
         % deathsChik
         deathsChik = pop(16);
        % deathsDen
         deathsDen = pop(17);
        %cumulativeInfectionChik
         incChik = pop(18);
        %cumulativeInfectionDen
         incDen = pop(19);
     
        dPop=zeros(19,1);

        %Total Population
         % %Total Human Population
        TotalHumanPopChik = sum(pop(1:5));
        TotalHumanPopDen = sum(pop(6:9));
        MtotChik = sum(pop(10:11))+1;
        MtotDen = sum(pop(12:13))+1;

        %Transmission variables
        Tchik= TranChik;
        gammachik = recChik;
        Tden= TranDen;
        gammaden = recDen;
        HRden = HRDen;
        lden=limRepeatInfDen;
        
        %Competing Mortality
        chikCompMort = compMortOn*mu(1)*(HRDen - 1)/TotalHumanPopChik; 
        denCompMort = compMortOn*mu(1)*(HRChik - 1)/TotalHumanPopDen; 


        %% Combined Model
        
        %Humans
        %S_Hchik S_Hden
        infRateChik = r*(Tchik(1,2)*I_Mchik)/TotalHumanPopChik;
        infRateDen = r*(Tden(1,2)*I_Mden)/TotalHumanPopDen;
        infRateDen2 = r*(Tden(1,2)*I_Mden)*lden/TotalHumanPopDen;
        recoverRateChik = gammachik(1);
        recoverRateDen = gammaden(1);
  
        %% Chikungunya model--- Susceptible Humans 0 Dengue
        %dS_Hchik/dt 
        dPop(1) = nu(1)*TotalHumanPopChik - infRateChik*S_Hchik - mu(1)*S_Hchik - chikCompMort*S_Hchik;  
        %I_H0chik
        dPop(2) = infRateChik*(1-symptChik)*(S_Hchik) - recoverRateChik*I_H0chik - mu(1)*I_H0chik - chikCompMort*I_H0chik; 
        %dI_H1chik/dt
        dPop(3)= infRateChik*symptChik*(S_Hchik)- recoverRateChik*I_H1chik - mu(1)*HRChik*I_H1chik - chikCompMort*I_H1chik;  
        %dR_H0chik/dt
        dPop(4)= recoverRateChik*(I_H0chik+(1-pSeqChik)*I_H1chik)+ pSeqRecovery*R_H1chik -mu(1)*R_H0chik - chikCompMort*R_H0chik; 
        %dR_H1chik/dt
        dPop(5)= recoverRateChik*pSeqChik*I_H1chik-pSeqRecovery*R_H1chik-mu(1)*R_H1chik - chikCompMort*R_H1chik; 
        
 
        %% Dengue Model
        %S_H0den/dt 
        dPop(6) = nu(1)*TotalHumanPopDen - infRateDen*S_H0den - mu(1)*S_H0den - denCompMort*S_H0den; 
        %S_H1den/dt
        dPop(7) =recoverRateDen*(I_H0den + I_H1den) - infRateDen2*S_H1den - mu(1)*S_H1den- denCompMort*S_H1den;   
        %I_H0den/dt
        dPop(8) =  infRateDen*(1-symptDen0)*S_H0den+ infRateDen2*S_H1den*(1-symptDen1) -recoverRateDen*I_H0den - mu(1)*I_H0den- denCompMort*I_H0den;  
        %I_H1den/dt
        dPop(9) =  infRateDen*S_H0den*symptDen0 + infRateDen2*S_H1den*symptDen1  -recoverRateDen*I_H1den - mu(1)*HRDen*I_H1den- denCompMort*I_H1den; 
        
        
        %% Mosquitos
        infChikMos = r*(Tchik(2,1)*(I_H0chik +I_H1chik))/TotalMosPopChik;
        infDenMos = r*(Tden(2,2)*I_Mden + Tden(2,1)*(I_H0den +I_H1den))/TotalMosPopDen;
        
        %% Chikungunya Transmission Mosquitos
        %dS_Mchik/dt
        dPop(10)=nu(2)*TotalMosPopChik -infChikMos*S_Mchik - mu(2)*S_Mchik;
        %dI_Mchik/dt
        dPop(11)= infChikMos*S_Mchik- mu(2)*I_Mchik;
        
        %% Dengue Transmission Mosquitos
        %dS_Mden/dt
        dPop(12)=nu(2)*TotalMosPopDen - infDenMos*S_Mden - mu(2)*S_Mden;
        %I_Mden/dt
        dPop(13)= infDenMos*S_Mden- mu(2)*I_Mden;

        %% Extra data
        %deathsTotChikModel
         dPop(14) = mu(1)*(S_Hchik + I_H0chik + HRChik*I_H1chik + R_H0chik + R_H1chik);
        %deathsTotDenModel
         dPop(15) = mu(1)*(S_H0den + S_H1den + I_H0den + HRDen*I_H1den);
         % deathsChik
         dPop(16) = mu(1)*(HRChik-1)*I_H1chik;
        % deathsDen
         dPop(17) = mu(1)*(HRDen-1)*I_H1den;
        %cumulativeInfectionChik
         dPop(18) = infRateChik*S_Hchik*symptChik;
        %cumulativeInfectionDen
         dPop(19) = infRateDen*S_H0den*symptDen0+infRateDen2*S_H1den*symptDen1;
      
    
    end


end
