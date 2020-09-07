function [t, population ] = combined_chikungunya_dengue_model_w_parallel_Zika(r,TranChik,TranDen,recChik,recDen,nu,mu, MaxTime, pSeqChik, pSeqRecovery, symptChik, symptDen0, symptDen1, HRChik, HRDen,limRepeatInfDen, vaxDenCoverage, vaxDenEfficacy, vaxChikEfficacy, denSens, denSpec, TranZika, recZika,symptZika, moreZikaMos, moreZikaPeople, zikaStart)
              
%Chikungunya Dengue Transmission models
%Anneke Claypool
%Version 2
%June 27, 2020

% This code runs the differential equations for the combined chikungunya
% and dengue models with a parallel Zika model, incorporating competing mortality. 

% This version includes the serotype test for dengue as a pre-rec for
% receiving the dengue vaccine. It also includes the sensitivity and
% specificity for that test. 

% Vaccination strategies are included, only vaccinate people who have not
% already been vaccinated. 
 
% This model is adapted from the MATLAB version 4.2 from page 123 of 
% "Modeling Infectious Disease in humans and animals" 
% by Keeling & Rohani.


%% Combined Model
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nsteps= MaxTime;
tspan = linspace(1,MaxTime, nsteps); %create timespan with results for each time step
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
compMortOn = 1;
VaxBounceBack = 1-0.6;

%% Initialize States

[t, population]=ode23t(@Diff,tspan,initialStates,options); % Use ode23t, The main ODE


%% Functions

% Calculates the differential rates used in the integration for chikungunya, dengue, and Zika.
    function dPop=Diff(t, pop)    
       
 % Susceptible Humans 0 Dengue
        S_HchikS_H0den = pop(1);
        I_H0chikS_H0den = pop(2);
        I_H1chikS_H0den = pop(3);
        R_H0chikS_H0den = pop(4);
        R_H1chikS_H0den = pop(5);

% Susceptible Humans 1 Dengue
        S_HchikS_H1den= pop(6);
        I_H0chikS_H1den= pop(7);
        I_H1chikS_H1den= pop(8);
        R_H0chikS_H1den = pop(9);
        R_H1chikS_H1den = pop(10);

% Infected Humans 0 Dengue
        S_HchikI_H0den = pop(11);
        I_H0chikI_H0den = pop(12);
        I_H1chikI_H0den = pop(13);
        R_H0chikI_H0den = pop(14);
        R_H1chikI_H0den = pop(15);

% Infected Humans 1 Dengue
        S_HchikI_H1den = pop(16);
        I_H0chikI_H1den = pop(17);
        I_H1chikI_H1den = pop(18);
        R_H0chikI_H1den = pop(19);
        R_H1chikI_H1den = pop(20);

%Total Human Population
        TotalHumanPop = sum(pop(1:20))+ sum(pop(30:39));

%Mosquitos
        S_MchikS_Mden = pop(21);
        I_MchikS_Mden = pop(22);
        S_MchikI_Mden = pop(23);
        I_MchikI_Mden = pop(24);
        
%Total Mosquito Population
        TotalMosPop = sum(pop(21:24));
        
% Extra states
        deathsTot = pop(25);
        deathsChik = pop(26);
        deathsDen = pop(27);
        incChik = pop(28);
        incDen = pop(29);
        
%Vaccine States
        %Dengue Vaccine
         
        S_HchikV_Hden = pop(30);
        I_H0chikV_Hden = pop(31);
        I_H1chikV_Hden = pop(32);
        R_H0chikV_Hden = pop (33);
        R_H1chikV_Hden = pop(34);        
        
        %Chikungunya Vaccine
        V_HchikS_H0den = pop(35);
        V_HchikS_H1den = pop(36);
        V_HchikI_H0den = pop(37);
        V_HchikI_H1den = pop(38);
        %Chikungunya and Dengue vaccine
        V_HchikV_Hden = pop(39);
        
        %% Zika!!
        %Add in infectious mosquitos at t = zikaStart
        addZikaInfMos = 0;
        addZikaPeople = 0;
            if (round(t,0) >= zikaStart && round(t,0) <= (zikaStart + 1)) 
                addZikaInfMos = moreZikaMos;
                addZikaPeople = moreZikaPeople; 
   
            end 
            
        %dS_Hzika/dt 
        S_Hzika = pop(42); 
        %I_H0zika
        I_H0zika = pop(43)+(addZikaPeople); 
        %dI_H1zika/dt
        I_H1zika = pop(44)+(addZikaPeople);  
        %dR_H0zika/dt
        R_Hzika = pop(45);
        %deathsTotZikaModel
        deathsTotZika = pop(46); 
        %cumulativeInfectionZika
        incZika = pop(47);

        % Zika Transmission Mosquito       
        %dS_Mzika/dt
        S_Mzika = pop(48);
        %I_Mzika/dt
        I_Mzika = pop(49)+addZikaInfMos;
        
        
%Total Human Population
        TotalHumanPop = sum(pop(1:20))+ sum(pop(30:39));

        %Chik initial variables
        S_Hchik = S_HchikS_H0den + S_HchikS_H1den + S_HchikI_H0den + S_HchikI_H1den + S_HchikV_Hden;
        S_Mchik = S_MchikS_Mden + S_MchikI_Mden;
        I_H0chik = I_H0chikS_H0den +  I_H0chikS_H1den + I_H0chikI_H0den+ I_H0chikI_H1den + I_H0chikV_Hden;
        I_H1chik = I_H1chikS_H0den + I_H1chikS_H1den + I_H1chikI_H0den +I_H1chikI_H1den + I_H1chikV_Hden; 
        I_Mchik= I_MchikS_Mden + I_MchikI_Mden;
        R_H0chik = R_H0chikS_H0den + R_H0chikS_H1den +R_H0chikI_H0den + R_H0chikI_H1den + R_H0chikV_Hden; 
        R_H1chik = R_H1chikS_H0den + R_H1chikS_H1den + R_H1chikI_H0den + R_H1chikI_H1den + R_H1chikV_Hden;
        
        %Den initial variables
        S_H0den = S_HchikS_H0den + I_H0chikS_H0den + I_H1chikS_H0den + R_H0chikS_H0den + R_H1chikS_H0den + V_HchikS_H0den;
        S_H1den = S_HchikS_H1den + I_H0chikS_H1den + I_H1chikS_H1den + R_H0chikS_H1den + R_H1chikS_H1den + V_HchikS_H1den; 
        S_Mden = S_MchikS_Mden + I_MchikS_Mden;
        I_H0den = S_HchikI_H0den + I_H0chikI_H0den + I_H1chikI_H0den +R_H0chikI_H0den + R_H1chikI_H0den + V_HchikI_H0den;
        I_H1den = S_HchikI_H1den + I_H0chikI_H1den + I_H1chikI_H1den + R_H0chikI_H1den +R_H1chikI_H1den + V_HchikI_H1den; 
        I_Mden= S_MchikI_Mden + I_MchikI_Mden;
        
        V_Hden = S_HchikV_Hden + I_H0chikV_Hden + I_H1chikV_Hden + R_H0chikV_Hden + R_H1chikV_Hden;
        V_Hchik = V_HchikS_H0den + V_HchikS_H1den + V_HchikI_H0den + V_HchikI_H1den;

        V_Hboth = V_HchikV_Hden;
        dPop=zeros(49,1);
        
        
        %Total Population
        Htot = sum(pop(1:20))+ sum(pop(30:39));
        Mtot = sum(pop(21:24));
        TotalMosPopZika = sum(pop(48:49));

        %Transmission variables
        Tchik= TranChik;
        gammachik = recChik;
        Tden= TranDen;
        gammaden = recDen;
        HRden = HRDen;
        lden=limRepeatInfDen;
        
        TotalHumanPopZika = sum(pop(42:45));
        Tzika= TranZika;
        gammazika = recZika;
        zikaCompMort = compMortOn*mu(1)*((HRChik - 1)+(HRDen - 1))/TotalHumanPopZika;


        %% Combined Model
        
        %Humans
        %S_Hchik S_Hden
        infRateChik = r*(Tchik(1,2)*I_Mchik)/Htot;
        infRateDen = r*(Tden(1,2)*I_Mden)/Htot;
        infRateDen2 = r*(Tden(1,2)*I_Mden)*lden/Htot;
        infRateZika = r*(Tzika(1,2)*I_Mzika)/TotalHumanPopZika;
        recoverRateChik = gammachik(1);
        recoverRateDen = gammaden(1);
        recoverRateZika = gammazika(1);
        
        %Vaccination Variables
        vaxRateWeek = vaxDenCoverage/52;
        sensitivity = denSens;
        specificity = denSpec;
        VaxProp = (vaxRateWeek*Htot)/(S_H0den + S_H1den);
        VaxPropS0 = VaxProp*(1-specificity);
        VaxPropS1 = VaxProp*(sensitivity);
        
        %% Susceptible Humans 0 Dengue
        %dS_Hchik S_H0den/dt 
        dPop(1) = nu(1)*TotalHumanPop - (infRateChik + infRateDen - infRateChik*infRateDen)*S_HchikS_H0den- VaxPropS0*S_HchikS_H0den - mu(1)*S_HchikS_H0den;  %checked
        %I_H0chik S_H0den
        dPop(2) = infRateChik*(1-symptChik)*(1-infRateDen)*(S_HchikS_H0den + (1-vaxChikEfficacy)*V_HchikS_H0den) - (infRateDen + recoverRateChik - infRateDen*recoverRateChik)*I_H0chikS_H0den - VaxPropS0*I_H0chikS_H0den- mu(1)*I_H0chikS_H0den; 
        %dI_H1chikS_H0den/dt
        dPop(3)= infRateChik*symptChik*(1-infRateDen)*(S_HchikS_H0den + (1-vaxChikEfficacy)*V_HchikS_H0den)- (infRateDen + recoverRateChik - infRateDen*recoverRateChik)*I_H1chikS_H0den - VaxPropS0*I_H1chikS_H0den- mu(1)*HRChik*I_H1chikS_H0den;  
        %dR_H0chikS_H0den/dt
        dPop(4)= recoverRateChik*(I_H0chikS_H0den+(1-pSeqChik)*I_H1chikS_H0den)*(1-infRateDen)+ pSeqRecovery*R_H1chikS_H0den*(1-infRateDen)- infRateDen*R_H0chikS_H0den - VaxPropS0*R_H0chikS_H0den-mu(1)*R_H0chikS_H0den;
        %dR_H1chikS_H0den/dt
        dPop(5)= recoverRateChik*pSeqChik*I_H1chikS_H0den*(1-infRateDen)-(pSeqRecovery + infRateDen - pSeqRecovery*infRateDen)*R_H1chikS_H0den- VaxPropS0*R_H1chikS_H0den-mu(1)*R_H1chikS_H0den;
        
        %% Susceptible Humans 1 Dengue
        %dS_HchikS_H1den/dt
        dPop(6) =recoverRateDen*(S_HchikI_H0den + S_HchikI_H1den)*(1-infRateChik) - (infRateChik + infRateDen2 - infRateChik*infRateDen2)*S_HchikS_H1den - VaxPropS1*S_HchikS_H1den - mu(1)*S_HchikS_H1den;  
        %dI_H0chikS_H1den/dt
        dPop(7) =recoverRateDen*(1-recoverRateChik)*(I_H0chikI_H0den + I_H0chikI_H1den) + infRateChik*(1-symptChik)*recoverRateDen*(S_HchikI_H0den + S_HchikI_H1den + (1-vaxChikEfficacy)*(V_HchikI_H0den+V_HchikI_H1den)) + infRateChik*(1-symptChik)*(S_HchikS_H1den+ (1-vaxChikEfficacy)*V_HchikS_H1den)*(1-infRateDen2) - (infRateDen2 +recoverRateChik - infRateDen2*recoverRateChik)*I_H0chikS_H1den - VaxPropS1*I_H0chikS_H1den - mu(1)*I_H0chikS_H1den; 
        %dI_H1chikS_H1den/dt
        dPop(8)= recoverRateDen*(1-recoverRateChik)*(I_H1chikI_H0den + I_H1chikI_H1den) + infRateChik*symptChik*recoverRateDen*(S_HchikI_H0den + S_HchikI_H1den+ (1-vaxChikEfficacy)*(V_HchikI_H0den+V_HchikI_H1den)) + infRateChik*symptChik*(S_HchikS_H1den+ (1-vaxChikEfficacy)*V_HchikS_H1den)*(1-infRateDen2) - (infRateDen2 + recoverRateChik -infRateDen2*recoverRateChik)*I_H1chikS_H1den - VaxPropS1*I_H1chikS_H1den - mu(1)*HRChik*I_H1chikS_H1den;  
        %dR_H0chikS_H1den/dt
        dPop(9)= recoverRateDen*(R_H0chikI_H0den + R_H0chikI_H1den) + recoverRateChik*(1-infRateDen2)*(I_H0chikS_H1den+(1-pSeqChik)*I_H1chikS_H1den)+ pSeqRecovery*R_H1chikS_H1den*(1-infRateDen2) + pSeqRecovery*recoverRateDen*(R_H1chikI_H0den + R_H1chikI_H1den) + recoverRateChik*recoverRateDen*(I_H0chikI_H0den + I_H0chikI_H1den + (1-pSeqChik)*(I_H1chikI_H0den+ I_H1chikI_H1den)) - infRateDen2*R_H0chikS_H1den - VaxPropS1*R_H0chikS_H1den -mu(1)*R_H0chikS_H1den;
        %dR_H1chikS_H1den/dt
        dPop(10)= recoverRateDen*(1-pSeqRecovery)*(R_H1chikI_H0den + R_H1chikI_H1den) + recoverRateChik*(1-infRateDen2)*pSeqChik*I_H1chikS_H1den + recoverRateChik*pSeqChik*recoverRateDen*(I_H1chikI_H0den + I_H1chikI_H1den) -(pSeqRecovery + infRateDen2 - pSeqRecovery*infRateDen2)*R_H1chikS_H1den- VaxPropS1*R_H1chikS_H1den -mu(1)*R_H1chikS_H1den;
       
        
        %% Infected Humans 0 Dengue        
        %dS_HchikI_H0den/dt
        dPop(11) =  (1-infRateChik)*infRateDen*(1-symptDen0)*S_HchikS_H0den+ (1-infRateChik)*infRateDen2*S_HchikS_H1den*(1-symptDen1) + (VaxBounceBack)*(1-infRateChik)*infRateDen2*(1-symptDen1)*(1-vaxDenEfficacy)*S_HchikV_Hden -(recoverRateDen + infRateChik - recoverRateDen*infRateChik)*S_HchikI_H0den - mu(1)*S_HchikI_H0den;  
        %dI_H0chikI_H0den/dt
        dPop(12) = infRateChik*(1-symptChik)*(S_HchikI_H0den + (1-vaxChikEfficacy)*V_HchikI_H0den)*(1-recoverRateDen) + (1-recoverRateChik)*infRateDen*I_H0chikS_H0den*(1-symptDen0)+ (1-recoverRateChik)*infRateDen2*I_H0chikS_H1den*(1-symptDen1)+ (1-recoverRateChik)*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*I_H0chikV_Hden*(1-symptDen1) + infRateChik*(1-symptChik)*infRateDen*(1-symptDen0)*(S_HchikS_H0den + (1-vaxChikEfficacy)*V_HchikS_H0den) +infRateChik*(1-symptChik)*(VaxBounceBack)*infRateDen2*(1-symptDen1)*((1-vaxDenEfficacy)*(1-vaxChikEfficacy)*V_HchikV_Hden+ (1-vaxDenEfficacy)*S_HchikV_Hden) + infRateChik*(1-symptChik)*infRateDen2*(1-symptDen1)*(S_HchikS_H1den+(1-vaxChikEfficacy)*V_HchikS_H1den) -(recoverRateDen + recoverRateChik - recoverRateDen*recoverRateChik)*I_H0chikI_H0den  - mu(1)*I_H0chikI_H0den; 
        %dI_H1chikI_H0den/dt
        dPop(13)=  infRateChik*symptChik*(S_HchikI_H0den+ (1-vaxChikEfficacy)*V_HchikI_H0den)*(1-recoverRateDen) + infRateChik*symptChik*infRateDen*(S_HchikS_H0den+ (1-vaxChikEfficacy)*V_HchikS_H0den)*(1-symptDen0)+ infRateChik*symptChik*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*S_HchikV_Hden*(1-symptDen1)+ infRateChik*symptChik*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*(1-vaxChikEfficacy)*V_HchikV_Hden*(1-symptDen1)+ infRateChik*symptChik*infRateDen2*(S_HchikS_H1den+ (1-vaxChikEfficacy)*V_HchikS_H1den)*(1-symptDen1) +(1-recoverRateChik)*infRateDen*I_H1chikS_H0den*(1-symptDen0)+(1-recoverRateChik)*(infRateDen2*(1-vaxDenEfficacy)*VaxBounceBack*I_H1chikV_Hden*(1-symptDen1) + infRateDen2*I_H1chikS_H1den*(1-symptDen1)) -(recoverRateDen + recoverRateChik - recoverRateDen*recoverRateChik)*I_H1chikI_H0den - mu(1)*HRChik*I_H1chikI_H0den;  
        %dR_H0chikI_H0den/dt
        dPop(14)=  recoverRateChik*(1-recoverRateDen)*(I_H0chikI_H0den+(1-pSeqChik)*I_H1chikI_H0den) + recoverRateChik*infRateDen*(1-symptDen0)*(I_H0chikS_H0den +(1-pSeqChik)*I_H1chikS_H0den)+ recoverRateChik*(VaxBounceBack)*infRateDen2*(1-symptDen1)*((1-vaxDenEfficacy)*I_H0chikV_Hden +(1-pSeqChik)*(1-vaxDenEfficacy)*I_H1chikV_Hden) + recoverRateChik*infRateDen2*(1-symptDen1)*(I_H0chikS_H1den+(1-pSeqChik)*I_H1chikS_H1den)+ pSeqRecovery*R_H1chikI_H0den*(1-recoverRateDen)+ infRateDen*R_H0chikS_H0den*(1-symptDen0)+ infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*R_H0chikV_Hden*(1-symptDen1)+ infRateDen2*R_H0chikS_H1den*(1-symptDen1)+ infRateDen*(1-symptDen0)*R_H1chikS_H0den*pSeqRecovery+ (VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*R_H1chikV_Hden*(1-symptDen1)*pSeqRecovery+ infRateDen2*R_H1chikS_H1den*pSeqRecovery*(1-symptDen1) -recoverRateDen*(R_H0chikI_H0den) -mu(1)*R_H0chikI_H0den;
        %dR_H1chikI_H0den/dt
        dPop(15)=  recoverRateChik*(pSeqChik*I_H1chikI_H0den)*(1-recoverRateDen) +(1-pSeqRecovery)*(infRateDen*R_H1chikS_H0den*(1-symptDen0)+ infRateDen2*R_H1chikS_H1den*(1-symptDen1)) +(1-pSeqRecovery)*(infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*R_H1chikV_Hden*(1-symptDen1)) + recoverRateChik*pSeqChik*(infRateDen*(1-symptDen0)*I_H1chikS_H0den + infRateDen2*(1-symptDen1)*I_H1chikS_H1den) + recoverRateChik*pSeqChik*infRateDen2*(VaxBounceBack)*(1-symptDen1)*(1-vaxDenEfficacy)*I_H1chikV_Hden-(recoverRateDen + pSeqRecovery - recoverRateDen*pSeqRecovery)*R_H1chikI_H0den - mu(1)*R_H1chikI_H0den;
        
        
        %% Infected Humans 1 Dengue
        %dS_HchikI_H1den/dt
        dPop(16) =  (1-infRateChik)*infRateDen*S_HchikS_H0den*symptDen0 + (1-infRateChik)*infRateDen2*S_HchikS_H1den*symptDen1 + (1-infRateChik)*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*S_HchikV_Hden*symptDen1 -(recoverRateDen + infRateChik -recoverRateDen*infRateChik)*S_HchikI_H1den  - mu(1)*HRDen*S_HchikI_H1den;  
        %dI_H0chikI_H1den/dt
        dPop(17) =  infRateChik*(1-symptChik)*(S_HchikI_H1den+ (1-vaxChikEfficacy)*V_HchikI_H1den)*(1-recoverRateDen) + infRateChik*(1-symptChik)*infRateDen*(S_HchikS_H0den + (1-vaxChikEfficacy)*V_HchikS_H0den)*symptDen0 + infRateChik*(1-symptChik)*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*S_HchikV_Hden*symptDen1 + infRateChik*(1-symptChik)*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*(1-vaxChikEfficacy)*V_HchikV_Hden*symptDen1 + infRateChik*(1-symptChik)*infRateDen2*(S_HchikS_H1den + (1-vaxChikEfficacy)*V_HchikS_H1den)*symptDen1 + (1- recoverRateChik)*infRateDen*I_H0chikS_H0den*symptDen0 + (1- recoverRateChik)*(infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*I_H0chikV_Hden*(symptDen1) + infRateDen2*I_H0chikS_H1den*symptDen1) -(recoverRateDen + recoverRateChik - recoverRateDen*recoverRateChik)*I_H0chikI_H1den - mu(1)*HRDen*I_H0chikI_H1den; 
        %dI_H1chikI_H1den/dt
        dPop(18)=   infRateChik*symptChik*S_HchikI_H1den*(1-recoverRateDen) + infRateChik*symptChik*(1-vaxChikEfficacy)*V_HchikI_H1den*(1-recoverRateDen) + infRateChik*(symptChik)*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*S_HchikV_Hden*symptDen1 + infRateChik*symptChik*(infRateDen*S_HchikS_H0den*symptDen0  + infRateDen*(1-vaxChikEfficacy)*V_HchikS_H0den*symptDen0+ infRateDen2*S_HchikS_H1den*symptDen1 + infRateDen2*(1-vaxChikEfficacy)*V_HchikS_H1den*symptDen1 + infRateDen2*(VaxBounceBack)*(1-vaxChikEfficacy)*(1-vaxDenEfficacy)*V_HchikV_Hden*symptDen1) + (1- recoverRateChik)*(infRateDen*I_H1chikS_H0den*symptDen0 + infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*I_H1chikV_Hden*symptDen1 + infRateDen2*I_H1chikS_H1den*symptDen1)  -(recoverRateDen + recoverRateChik- recoverRateDen*recoverRateChik)*I_H1chikI_H1den - mu(1)*HRDen*HRChik*I_H1chikI_H1den;  
        %dR_H0chikI_H1den/dt
        dPop(19)=   recoverRateChik*(1-recoverRateDen)*(I_H0chikI_H1den+(1-pSeqChik)*I_H1chikI_H1den) + infRateDen*R_H0chikS_H0den*symptDen0 + infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*R_H0chikV_Hden*symptDen1 + infRateDen2*R_H0chikS_H1den*symptDen1 + pSeqRecovery*R_H1chikI_H1den*(1-recoverRateDen) + pSeqRecovery*(infRateDen*symptDen0*R_H1chikS_H0den + (VaxBounceBack)*infRateDen2*symptDen1*(1-vaxDenEfficacy)*R_H1chikV_Hden +infRateDen2*symptDen1*R_H1chikS_H1den)   + recoverRateChik*(symptDen0*infRateDen*(I_H0chikS_H0den+(1-pSeqChik)*I_H1chikS_H0den)+ symptDen1*(VaxBounceBack)*infRateDen2*((1-vaxDenEfficacy)*I_H0chikV_Hden+(1-pSeqChik)*(1-vaxDenEfficacy)*I_H1chikV_Hden)  + infRateDen2*symptDen1*(I_H0chikS_H1den+(1-pSeqChik)*I_H1chikS_H1den))   -recoverRateDen*(R_H0chikI_H1den) -mu(1)*HRDen*R_H0chikI_H1den;
        %dR_H1chikI_H1den/dt
        dPop(20)=   recoverRateChik*pSeqChik*I_H1chikI_H1den*(1-recoverRateDen) + (1-pSeqRecovery)*(infRateDen*R_H1chikS_H0den*symptDen0 + infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*R_H1chikV_Hden*symptDen1 + infRateDen2*R_H1chikS_H1den*symptDen1) + recoverRateChik*pSeqChik*(infRateDen*I_H1chikS_H0den*symptDen0 + infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*I_H1chikV_Hden*symptDen1 + infRateDen2*I_H1chikS_H1den*symptDen1) -(recoverRateDen + pSeqRecovery - recoverRateDen*pSeqRecovery)*R_H1chikI_H1den -mu(1)*HRDen*R_H1chikI_H1den;
        
        %% Vaccinated Humans Dengue
        %dS_HchikV_Hden/dt
        dPop(30) = VaxPropS0*S_HchikS_H0den + VaxPropS1*S_HchikS_H1den -(VaxBounceBack*infRateDen2*(1-vaxDenEfficacy) + infRateChik - VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*infRateChik)*S_HchikV_Hden  - mu(1)*S_HchikV_Hden;  
        %dI_H0chikV_Hden/dt
        dPop(31) = VaxPropS0*I_H0chikS_H0den + VaxPropS1*I_H0chikS_H1den + (1-vaxChikEfficacy)*infRateChik*(1-symptChik)*(1-(infRateDen2*(1-vaxDenEfficacy)))*V_HchikV_Hden + infRateChik*(1-symptChik)*(1-(infRateDen2*(1-vaxDenEfficacy)))*S_HchikV_Hden-(VaxBounceBack*infRateDen2*(1-vaxDenEfficacy) + recoverRateChik - VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*recoverRateChik)*I_H0chikV_Hden - mu(1)*I_H0chikV_Hden; 
        %%dI_H1chikV_Hden/dt
        dPop(32) = VaxPropS0*I_H1chikS_H0den + VaxPropS1*I_H1chikS_H1den + (1-vaxChikEfficacy)*infRateChik*(symptChik)*(1-(infRateDen2*(1-vaxDenEfficacy)))*V_HchikV_Hden +infRateChik*symptChik*(1-(infRateDen2*(1-vaxDenEfficacy)))*S_HchikV_Hden-(VaxBounceBack*infRateDen2*(1-vaxDenEfficacy) + recoverRateChik - VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*recoverRateChik)*I_H1chikV_Hden- mu(1)*HRChik*I_H1chikV_Hden;  
        %dR_H0chikV_Hden/dt
        dPop(33) = VaxPropS0*R_H0chikS_H0den + VaxPropS1*R_H0chikS_H1den + recoverRateChik*(1-(infRateDen2*(1-vaxDenEfficacy)))*I_H0chikV_Hden + recoverRateChik*(1-pSeqChik)*(1-(infRateDen2*(1-vaxDenEfficacy)))*I_H1chikV_Hden + pSeqRecovery*(1-(infRateDen2*(1-vaxDenEfficacy)))*R_H1chikV_Hden -VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*R_H0chikV_Hden -mu(1)*R_H0chikV_Hden;
        %dR_H1chikV_Hden/dt
        dPop(34) = VaxPropS0*R_H1chikS_H0den + VaxPropS1*R_H1chikS_H1den + recoverRateChik*pSeqChik*(1-(infRateDen2*(1-vaxDenEfficacy)))*I_H1chikV_Hden-(VaxBounceBack*infRateDen2*(1-vaxDenEfficacy) + pSeqRecovery - VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*pSeqRecovery)*R_H1chikV_Hden -mu(1)*R_H1chikV_Hden;
        
        %% Vaccinated Humans Chikungunya
        %dV_Hchik S_H0den/dt 
        dPop(35) =  -((1-vaxChikEfficacy)*infRateChik + infRateDen - (1-vaxChikEfficacy)*infRateChik*infRateDen)*V_HchikS_H0den -VaxPropS0*V_HchikS_H0den - mu(1)*V_HchikS_H0den;  
        %dV_HchikS_H1den/dt
        dPop(36) = (1-(infRateChik*(1-vaxChikEfficacy)))*recoverRateDen*V_HchikI_H0den + (1-(infRateChik*(1-vaxChikEfficacy)))*recoverRateDen*V_HchikI_H1den- ((1-vaxChikEfficacy)*infRateChik + infRateDen2 - ((1-vaxChikEfficacy)*infRateChik*infRateDen2))*V_HchikS_H1den -VaxPropS1*V_HchikS_H1den - mu(1)*V_HchikS_H1den;  
        %dV_HchikI_H0den/dt
        dPop(37) =  (1-(infRateChik*(1-vaxChikEfficacy)))*VaxBounceBack*infRateDen2*(1-symptDen1)*(1-vaxDenEfficacy)*V_HchikV_Hden + (1-(infRateChik*(1-vaxChikEfficacy)))*infRateDen*(1-symptDen0)*V_HchikS_H0den + (1-(infRateChik*(1-vaxChikEfficacy)))*infRateDen2*(1-symptDen1)*V_HchikS_H1den-(recoverRateDen + (1-vaxChikEfficacy)*infRateChik - (recoverRateDen*(1-vaxChikEfficacy)*infRateChik))*V_HchikI_H0den - mu(1)*V_HchikI_H0den;  
        %dV_HchikI_H1den/dt
        dPop(38) =  (1-(infRateChik*(1-vaxChikEfficacy)))*VaxBounceBack*infRateDen2*(symptDen1)*(1-vaxDenEfficacy)*V_HchikV_Hden + (1-(infRateChik*(1-vaxChikEfficacy)))*infRateDen*symptDen0*V_HchikS_H0den + (1-(infRateChik*(1-vaxChikEfficacy)))*infRateDen2*symptDen1*V_HchikS_H1den -(recoverRateDen + (1-vaxChikEfficacy)*infRateChik -(recoverRateDen*(1-vaxChikEfficacy)*infRateChik))*V_HchikI_H1den  - mu(1)*HRDen*V_HchikI_H1den;  
        %dV_HchikV_Hden/dt
        dPop(39) = VaxPropS0*V_HchikS_H0den + VaxPropS1*V_HchikS_H1den -(VaxBounceBack*infRateDen2*(1-vaxDenEfficacy) + (1-vaxChikEfficacy)*infRateChik -VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*(1-vaxChikEfficacy)*infRateChik)*V_HchikV_Hden  - mu(1)*V_HchikV_Hden;  

        
        %% Mosquitos
        infChikMos = r*(Tchik(2,1)*(I_H0chik +I_H1chik))/Mtot;
        infDenMos = r*(Tden(2,2)*I_Mden + Tden(2,1)*(I_H0den +I_H1den))/Mtot;
        infZikaMos = r*(Tzika(2,1)*(I_H0zika +I_H1zika))/TotalMosPopZika;

        %dS_MchikS_Mden/dt
        dPop(21)=nu(2)*TotalMosPop - (infChikMos + infDenMos -infChikMos*infDenMos)*S_MchikS_Mden - mu(2)*S_MchikS_Mden;
        %dI_MchikS_Mden/dt
        dPop(22)= infChikMos*(1-infDenMos)*S_MchikS_Mden-infDenMos*I_MchikS_Mden- mu(2)*I_MchikS_Mden;
        %dS_MchikI_Mden/dt
        dPop(23)= infDenMos*(1-infChikMos)*S_MchikS_Mden-infChikMos*S_MchikI_Mden- mu(2)*S_MchikI_Mden;
        %dI_MchikI_Mden/dt
        dPop(24)= infDenMos*I_MchikS_Mden + infChikMos*S_MchikI_Mden + infChikMos*infDenMos*S_MchikS_Mden- mu(2)*I_MchikI_Mden;
        
        
        %% Extra data
        %deathsTot
         dPop(25) = mu(1)*(S_HchikS_H0den + I_H0chikS_H0den  + I_H1chikS_H0den*HRChik  + R_H0chikS_H0den  + R_H1chikS_H0den  + S_HchikS_H1den + I_H0chikS_H1den + I_H1chikS_H1den*HRChik + R_H0chikS_H1den  + R_H1chikS_H1den  + S_HchikI_H0den  + I_H0chikI_H0den  + I_H1chikI_H0den*HRChik  + R_H0chikI_H0den + R_H1chikI_H0den  + HRDen*(S_HchikI_H1den  + I_H0chikI_H1den  + I_H1chikI_H1den*HRChik  + R_H0chikI_H1den  + R_H1chikI_H1den)+S_HchikV_Hden + I_H0chikV_Hden + HRChik*I_H1chikV_Hden + R_H0chikV_Hden + R_H1chikV_Hden +V_HchikS_H0den +  V_HchikS_H1den + V_HchikI_H0den + HRDen*V_HchikI_H1den + V_HchikV_Hden);
        % deathsChik
         dPop(26) = mu(1)*(HRChik-1)*(I_H1chikS_H0den  + I_H1chikS_H1den + I_H1chikI_H0den + I_H1chikI_H1den + I_H1chikV_Hden);
        % deathsDen
         dPop(27) = mu(1)*(HRDen-1)*(S_HchikI_H1den + I_H0chikI_H1den  + I_H1chikI_H1den  + R_H0chikI_H1den  + R_H1chikI_H1den + V_HchikI_H1den);
        %cumulativeInfectionChik
         dPop(28) = (infRateChik*S_Hchik*symptChik)+ V_Hchik*infRateChik*symptChik*(1-vaxChikEfficacy) + V_Hboth*infRateChik*symptChik*(1-vaxChikEfficacy);
        %cumulativeInfectionDen
         dPop(29) = infRateDen*S_H0den*symptDen0+infRateDen2*S_H1den*symptDen1+ V_Hden*infRateDen2*symptDen1*(1-vaxDenEfficacy)+ V_Hboth*infRateDen2*symptDen1*(1-vaxDenEfficacy);
      
        % Number of people tested
        totalTests = (vaxRateWeek*Htot);
        dPop(40)=totalTests;
        % Number of people vaccinated
        totalVax = VaxPropS0*S_HchikS_H0den + VaxPropS1*S_HchikS_H1den +VaxPropS0*I_H0chikS_H0den + VaxPropS1*I_H0chikS_H1den + VaxPropS0*I_H1chikS_H0den + VaxPropS1*I_H1chikS_H1den +VaxPropS0*R_H0chikS_H0den + VaxPropS1*R_H0chikS_H1den +VaxPropS0*R_H1chikS_H0den + VaxPropS1*R_H1chikS_H1den +  VaxPropS0*V_HchikS_H0den + VaxPropS1*V_HchikS_H1den;
        dPop(41)=totalVax;
        
        %% Zika model--- Susceptible Humans 0 Dengue
        %dS_Hzika/dt 
        dPop(42) = nu(1)*TotalHumanPopZika - infRateZika*S_Hzika - mu(1)*S_Hzika - zikaCompMort*S_Hzika; 
        %I_H0zika
        dPop(43) = infRateZika*(1-symptZika)*(S_Hzika) - recoverRateZika*I_H0zika - mu(1)*I_H0zika - zikaCompMort*I_H0zika; 
        %dI_H1zika/dt
        dPop(44)= infRateZika*symptZika*(S_Hzika)- recoverRateZika*I_H1zika - mu(1)*I_H1zika - zikaCompMort*I_H1zika; 
        %dR_Hzika/dt
        dPop(45)= recoverRateZika*(I_H0zika+I_H1zika)-mu(1)*R_Hzika - zikaCompMort*R_Hzika;

         %% Extra Zika data
        %deathsTotZikaModel
         dPop(46) = (mu(1)+zikaCompMort)*(S_Hzika + I_H0zika + I_H1zika + R_Hzika);
        %cumulativeInfectionZika
         dPop(47) = infRateZika*S_Hzika*symptZika;
        %% Zika Transmission Mosquitos
        %dS_Mzika/dt
        dPop(48)=nu(2)*TotalMosPopZika -infZikaMos*S_Mzika - mu(2)*S_Mzika; 
        %dI_Mzika/dt
        dPop(49)= addZikaInfMos + infZikaMos*S_Mzika- mu(2)*I_Mzika;
        

    
    end


end
