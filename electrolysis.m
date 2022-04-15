function [mdot6, mdot7, h5, s5, h6, s6, h7, s7, T5, T6, T7, p5, p6, p7, Erev, Jm, Wdote, X6, X7, Y6, Y7] = electrolysis(NM, M, mdot5, dhn, dsn, p1, T1)
%% Initialize REFPROP10
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
%% Set up system
mdot6 = zeros(1,NM);
mdot7 = zeros(1,NM);
Y6 = zeros(1, NM);
Y7 = zeros(1, NM);
h5m = zeros(1,NM);
s5m = zeros(1,NM);
g5m = zeros(1,NM);
h6m = zeros(1,NM);
s6m = zeros(1,NM);
g6m = zeros(1,NM);
h7m = zeros(1,NM);
s7m = zeros(1,NM);
g7m = zeros(1,NM);
n = 2;                      % Number of electrons transferred per formula conversion
F = 96485.3329;             % [A*s/mol] Faraday constant
p5 = p1;                    % Pa
T5 = T1;                    % K
p6 = p5;                    % Pa
T6 = T5;                    % K
p7 = p5;                    % Pa
T7 = T5;                    % K
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
% Flow rates
ndot5 = mdot5(8)/(M(1)/1000);  % Molar flow rate [mol/s]
% H2O => H2 + 0.5*O2
ndot6 = 0.5*ndot5;
mdot6(3) = ndot6*(M(3)/1000);  % Mass flow rate of oxygen [kg/s]
ndot7 = ndot5;
mdot7(2) = ndot7*(M(2)/1000);   % Mass flow rate of hydrogen [kg/s]

% Find the mass fractions at States 6 and 7
for i = 1:NM
    Y6(i) = mdot6(i)/sum(mdot6);
    Y7(i) = mdot7(i)/sum(mdot7);
end
% Find the mole fractions at States 6 and 7
X6 = Y2X(Y6, M);
X7 = Y2X(Y7, M);

%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
sm5 = RP.SETFLUIDSdll('WATER.FLD');
r5 = RP.ABFLSHdll('TP', T5, p5/1000, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
h5m(8) = r5.h + dhn(8);  % [J/mol]
s5m(8) = r5.s + dsn(8);  % [J/molK]
g5m(8) = h5m(8) - T5*s5m(8);    % [J/mol]
sm6 = RP.SETFLUIDSdll('OXYGEN.FLD');
r6 = RP.ABFLSHdll('TP', T6, p6/1000, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
h6m(3) = r6.h + dhn(3);  % [J/mol]
s6m(3) = r6.s + dsn(3);  % [J/molK]
g6m(3) = h6m(3) - T6*s6m(3);    % [J/mol]
sm7 = RP.SETFLUIDSdll('HYDROGEN.FLD');
r7 = RP.ABFLSHdll('TP', T7, p7/1000, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
h7m(2) = r7.h + dhn(2);  % [J/mol]
s7m(2) = r7.s + dsn(2);  % [J/molK]
g7m(2) = h7m(2) - T7*s7m(2);    % [J/mol]
% Calculate the Reversible Cell Potential
Erev = (g5m(8) - g7m(2) - 0.5*g6m(3))/(n*F);        % Volts
Jm = abs(Erev)*n*F;                                 % J/molH2O
Wdote = Jm*ndot5;                                 % Energy [W] required at the water flow rate specified
%% Convert enthalpies and entropies to J/kg
for i=1:NM
    h5(i) = h5m(i)/(M(i)/1000);                 % [J/mol] / ([gm/mol] * [kg/1000 gm]) = [J/kg]
    s5(i) = s5m(i)/(M(i)/1000);                 % [J/molK] / ([gm/mol] * [kg/1000 gm]) = [J/kgK]
    h6(i) = h6m(i)/(M(i)/1000);                 % [J/mol] / ([gm/mol] * [kg/1000 gm]) = [J/kg]
    s6(i) = s6m(i)/(M(i)/1000);                 % [J/molK] / ([gm/mol] * [kg/1000 gm]) = [J/kgK]    
    h7(i) = h7m(i)/(M(i)/1000);                 % [J/mol] / ([gm/mol] * [kg/1000 gm]) = [J/kg]
    s7(i) = s7m(i)/(M(i)/1000);                 % [J/molK] / ([gm/mol] * [kg/1000 gm]) = [J/kgK]
end