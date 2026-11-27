function [k10, k12, k13, k21, k31, k1e, ke0, Vc, gamma, E0, Emax, C50, Weight]=propofol_patient_parameters(Id,Patients)
%% Demographic Schnider PK parameters and the tabulated PD parameters.
% Input: profile index and the tabulated patient-profile matrix.
% Output: nominal Schnider PK rates/volume plus tabulated Hill parameters.
% Ce50, gamma, E0 and Emax are read directly from the profile table; they
% are not inferred from demographic data by the Schnider population model.

%% DATA

% Read the demographic columns of the supplied reference-profile matrix.
Age    = Patients( Id , 2 );
Height = Patients( Id , 3 );
Weight = Patients( Id , 4 );
Gender = Patients( Id , 5 );

%% MODEL PARAMETERS

if Gender == 1
    %lean body mass (James Formula for Men)
    lbm = 1.1*Weight - 128*( Weight / Height )^2;
else % Gender code 2: women
    %lean body mass (James Formula for Women)
    lbm = 1.07*Weight - 148*( Weight / Height )^2;
end

%Volume of the compartments (Schnider model for propofol)
Vc = 4.27;                                                          %[L]
V2 = 18.9 - 0.391*( Age - 53 );                                     %[L]
V3 = 238;                                                          %[L]

%Clearances are in L/min; division by 60 below converts rate constants to 1/s.
Cl1 = 1.89 + 0.0456*( Weight - 77 ) - 0.0681*( lbm - 59 ) +...      
      0.0264*( Height - 177 );                                      %[L/min]
Cl2 = 1.29 - 0.024*( Age - 53 );                                    %[L/min]
Cl3 = 0.836;                                                        %[L/min]

%Transfer Rate of Drug
%N.B. kij transfer rate from i to j
k10	= Cl1 / Vc / 60;                                                %[1/s]
k12	= Cl2 / Vc / 60;                                                %[1/s]
k13	= Cl3 / Vc / 60;                                                %[1/s]
k21	= Cl2 / V2 / 60;                                                %[1/s]
k31	= Cl3 / V3 / 60;                                                %[1/s]
k41	= 0.459 / 60;                                                   %[1/s]
ke0  = k41;                                                         %[1/s]
%equivalence hypothesis
k1e = ke0;


%% Hill function parameters

C50 = Patients( Id, 6 );
gamma = Patients( Id, 7 );
E0 = Patients( Id, 8);
Emax = Patients( Id, 9);

end
