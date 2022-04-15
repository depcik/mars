%%********************************************************************%%
% 
% Authors: Dr. Christopher Depcik (depcik@ku.edu)
%
% This is written as a long file, but it might be good to break it up into 
% subroutines that can be called
% 
% At each state, let's have values of: pressure [Pa], temperature [K], 
% and arrays of enthalpies [J/kg] and entropies [J/kgK]
% Moreover, we will need to adjust all of the states to use the same datum
% as per the NIST webbook
%%********************************************************************%%
%% Input parameters
clear
VwtrTank = 100*0.001;           % Size of water tank: L=>m3
VCO2Tank = 100*0.001;           % Size of CO2 tank: L=>m3
NM = 9;                         % Number of chemical species
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
Xatm = [0, 0, 0.00174, 0.0208, 0.9490, 0.000747, 0.0279, 0, 0];       % Atmosphere on Mars [by volume]
pMars = 610;                    % Martian atmospheric pressure [Pa]
TMars = 210;                    % Maritan atmospheric temperature [K]
%% Initialize REFPROP10
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
%% Common Variables
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
%   https://webbook.nist.gov/
M = [18.01528, 2*1.00784, 2*15.999, 39.948, 44.01, 28.01, 2*14.0067, 18.01528, 16.04];   % Molecular Masses [gm/mol]
% Calculate the normalization parameters
[dhn, dsn] = speciesnorm();     % Enthalpy returns [J/mol], entropy [J/molK]
%% Mars Water Mining
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
h1 = zeros(1, NM);
s1 = zeros(1, NM);
sm = RP.SETFLUIDSdll('WATER.FLD');
mdot1 = zeros(1,NM);
mdot1(8) = 200/40/1000/60;          % H2O(l) mass flow rate out of mining [kg/s]
Y1 = zeros(1, NM);                  % Species mass fractions at State 1
for i = 1:NM
    Y1(i) = mdot1(i)/sum(mdot1);
end
X1 = Y2X(Y1, M);                    % Mole fractions at State 1 [-]
T(1) = 273.15;                      % Temperature of water in tank [K]
p(1) = 0.13523*1e6;                 % Pressure of water in tank [Pa]
Wdot1 = 340;                        % Power it takes for the mass flow rate determined [W]
% Check the water in the tank to make sure it is a liquid
r1 = RP.ABFLSHdll('TP', T(1), p(1)/1000, zu, int8(0));
% r.h [J/mol]
% r.s [J/molK]
h1(8) = (r1.h + dhn(8))/M(8)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s1(8) = (r1.s + dsn(8))/M(8)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
if (r1.q == -998)
    'Water Tank - Subcooled Liquid'          % It should be a subcooled liquid for storage
end
%% Mars Atmosphere CO2 Collection
mdot2 = zeros(1,NM);
mdot2(5) = 20/1000/3600;        % Mass flow rate out of CO2 collection [kg/s]
eta2 = 0.50;                    % CO2 capture efficiency [-]
T(2) = 273.15;                    % Temperature of CO2 in tank [K]
p(2) = 3.4852*1e6;                % Pressure of CO2 in tank [Pa]
T(3) = TMars;                       % Temperature of Mars atmosphere [K]
p(3) = pMars;                       % Pressure of Mars atmosphere [Pa]
Wdot2 = 50;                     % Power it takes for the mass flow rate determined [W]
[h2, s2, mdot3, h3, s3, mdot4, h4, s4, p(4), T(4), X2, X3, X4, Y2, Y3, Y4] = atmCO2col(NM, M, Xatm, mdot2, T, p, dhn, dsn, eta2);
%% Electrolysis of Liquid Water
mdot5 = zeros(1,NM);
mdot5(8) = 1.1/1000;        % Flow rate out of the Water tank [kg/s] - the rest of the system is based off this
Y5 = zeros(1, NM);          % Mass fractions of gas species at State 5
for i = 1:NM
    Y5(i) = mdot5(i)/sum(mdot5);
end
X5 = Y2X(Y5, M); % Mole fractions of gas species at State 5
[mdot6, mdot7, h5, s5, h6, s6, h7, s7, T(5), T(6), T(7), p(5), p(6), p(7), Erev, Jm, Wdot567, X6, X7, Y6, Y7] = electrolysis(NM, M, mdot5, dhn, dsn, p(1), T(1));
%% O2 Liquefaction
% First thing is to find the minimum enthalpy at the temperature specifed
etac = 0.70;            % Compressor isentropic efficiency [-]
% TO2store = 203.15;
TO2store = 210; % Setting this to 210 K to match the manuscript values
[mdot8, mdot9, mdot10, mdot11, mdot12, mdot13, mdot14, p, T, h8, h9, h10, ...
    h11, h12, h13, h14, s8, s9, s10, s11, s12, s13, s14, q12, Wdot910, Qdot1011, ...
    Qdot1314, X8, X9, X10, X11, X12, X13, X14, Y8, Y9, Y10, Y11, Y12, Y13, Y14, sigmadot910, sigmadot1011, sigmadot1112, sigmadot1314, sigmadot6149] = O2Liquefaction(NM, M, TO2store, p, T, mdot6, h6, s6, etac, dhn, dsn);
%% Sabatier Process
% 4H2 + CO2 => CH4 + 2H2O
% Molar flow rate ratio: 4 to 1 (H2 to CO2)
H2Cm = 4;
etaCO2 = 0.58;          % CO2 conversion efficiency
SCH4 = 0.995;           % CH4 selection efficiency
pSab = 98066.5;         % Pressure for Sabatier process [Pa] 
% Note regarding pressure: this will depend on the Sabatier catalyst chosen
% and the reactor conditions. According to Le Chatelier, we want the
% pressure high to promote CH4 generation. Currently, we are going to
% ignore this pressure and let it be the same as that exiting the
% Electrolysis reaction, which is the same as the water tank. This might
% have to be revisited later
TSab = 583.15;          % Temperature for Sabatier process [K]
etac = 0.70;            % Compressor efficiency [-]
p(22) = 1.82385*1e6;    % Inlet pressure to hydrogen/CO2 separator [Pa]
T(22) = 294.15;         % Inlet temperature to hydrogen/CO2 separator [K]
etaH2rec = 0.871;       % Hydrogen recovery from separator on a mass basis
etaCO2rec = 0.940;      % CO2 recovery from separator on a mass basis
[mdot15, mdot16, mdot17, mdot18, mdot19, mdot20, mdot21, mdot22, ...
    mdot23, mdot24, mdot25, mdot26, mdot27, p, T, h15, h16, h17, h18, h19, h20, ...
    h21, h22, h23, h24, h25, h26, h27, s15, s16, s17, s18, s19, s20, s21, s22, s23, ...
    s24, s25, s26, s27, Qdot1516, Qdot717, Qdot161718, Qdot181920, Wdot2021, ...
    Qdot2122, Wdot22232425, Qdot72617, Wdot2427, Qdot152716,...
    Y15, Y16, Y17, Y18, Y19, Y20, Y21, Y22, Y23, Y24, Y25, Y26, Y27,...
    X15, X16, X17, X18, X19, X20, X21, X22, X23, X24, X25, X26, X27] = ...
    Sabatier(NM, M, H2Cm, etaCO2, SCH4, pSab, TSab, etac, etaH2rec, etaCO2rec, p, T, ...
    mdot7, h7, h2, s2, dhn, dsn);
%% Find the mass fractions leaving the CO2/H2 separator
% These are already found in the Sabatier block now.
% Based on the pressure of State 25, what would be the temperatures at
% which H2, CO2, CO, and CH4 liquefy?
% Example at 0.101325 MPa
% Find Sublimation Temperature of CO2
%sm = RP.SETFLUIDSdll('CO2.FLD');
%rtst = RP.ABFLSHdll('PQ', 101.325, -98, zu, int8(0));
%Tfrz = rtst.T;
%% Methane Liquefaction
% OK, let's first find the pressure that generates the lowest enthalpy of
% methane given the coldest temperature we can get it to using the Martian atmosphere
% So, now we need to find the pressure that provides the smallest value of
% enthalpy through the valve
% TCH4store = 203.15;
TCH4store = 210; % Setting this to 210 K to match the manuscript values
etac = 0.70;
%% Initialization
mdot28 = zeros(1,NM);
mdot29 = zeros(1,NM);
mdot30 = zeros(1,NM);
mdot31 = zeros(1,NM);
mdot32 = zeros(1,NM);
Y28 = zeros(1, NM);
Y29 = zeros(1, NM);
Y30 = zeros(1, NM);
Y31 = zeros(1, NM);
Y32 = zeros(1, NM);
h28 = zeros(1,NM);
h29 = zeros(1,NM);
h30 = zeros(1,NM);
h31 = zeros(1,NM);
h32 = zeros(1,NM);
s28 = zeros(1,NM);
s29 = zeros(1,NM);
s30 = zeros(1,NM);
s31 = zeros(1,NM);
s32 = zeros(1,NM);
%% Compress to Achieve Liquid CO2 for Tanks
% Conservation of Mass
mdot28 = mdot25;
% Adding Dalton's partial pressure functionality
% Mass fractions
for i = 1:NM
    Y28(i) = mdot28(i)/sum(mdot28);
end
% Mole fractions
X28 = Y2X(Y28, M);
% Conservation of Momentum
p(28) = p(2);
% Conservation of Energy
% Need to solve isentropic first and use the compressor efficiency
% Conservation of Entropy
m25s25 = 0;
m25h25 = 0;
for i=1:NM
    m25s25 = m25s25 + mdot25(i)*s25(i);     % [W/K]
    m25h25 = m25h25 + mdot25(i)*h25(i);     % [W]
end
T28g = T(25)*((p(28)/p(25))^(1.25-1));  % Guess of isentropic temperature: State 21s [K]
T28s = fminsearch(@objectivefcn2, T28g, [], p(28), m25s25, mdot28, M, dsn, NM, X28); % X28 for partial pressures
% Now, we need to find the enthalpy at State 21s; PARTIAL PRESSURES USED TO
% QUERY STATE PROPERTIES
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
h28s = zeros(1,NM);
sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
r28s = RP.ABFLSHdll('TP', T28s, X28(2)*p(28)/1000, zu, int8(0));
h28s(2) = (r28s.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
sm = RP.SETFLUIDSdll('CO2.FLD');
r28s = RP.ABFLSHdll('TP', T28s, X28(5)*p(28)/1000, zu, int8(0));
h28s(5) = (r28s.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
sm = RP.SETFLUIDSdll('METHANE.FLD');
r28s = RP.ABFLSHdll('TP', T28s, X28(9)*p(28)/1000, zu, int8(0));
h28s(9) = (r28s.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
sm = RP.SETFLUIDSdll('CO.FLD');
r28s = RP.ABFLSHdll('TP', T28s, X28(6)*p(28)/1000, zu, int8(0));
h28s(6) = (r28s.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
m28h28s = 0;
for i=1:NM
    m28h28s = m28h28s + mdot28(i)*h28s(i);     % [W]
end
% Use the isentropic efficiency equation
m28h28 = (m25h25*etac + (m28h28s-m25h25))/etac;   % [W]
% Now, we have to find a temperature that achieves this enthalpy
T(28) = fminsearch(@objectivefcn3, T28s, [], p(28), m28h28, mdot28, M, dhn, NM, X28); % Partial pressures used here
% Finally, we have State 21 defined and we can calculate the power needed;
% PARTIAL PRESSURES HAVE BEEN ADDED
sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
r28 = RP.ABFLSHdll('TP', T(28), X28(2)*p(28)/1000, zu, int8(0));
h28(2) = (r28.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s28(2) = (r28.s + dsn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
sm = RP.SETFLUIDSdll('CO2.FLD');
r28 = RP.ABFLSHdll('TP', T(28), X28(5)*p(28)/1000, zu, int8(0));
h28(5) = (r28.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s28(5) = (r28.s + dsn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
sm = RP.SETFLUIDSdll('METHANE.FLD');
r28 = RP.ABFLSHdll('TP', T(28), X28(9)*p(28)/1000, zu, int8(0));
h28(9) = (r28.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s28(9) = (r28.s + dsn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
sm = RP.SETFLUIDSdll('CO.FLD');
r28 = RP.ABFLSHdll('TP', T(28), X28(6)*p(28)/1000, zu, int8(0));
h28(6) = (r28.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s28(6) = (r28.s + dsn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
Wdot2528 = m25h25 - m28h28;
m28s28 = 0;
for i = 1:NM
    m28s28 = m28s28 + mdot28(i) * s28(i);
end
sigmadot2528 = m28s28 - m25s25;
%% Separate out Liquid CO2
% Conservation of Mass
mdot29 = zeros(1,NM);
mdot30 = zeros(1,NM);
h29 = zeros(1,NM);
s29 = zeros(1,NM);
h30 = zeros(1,NM);
s30 = zeros(1,NM);
h31 = zeros(1, NM);
h32 = zeros(1, NM);
h33 = zeros(1, NM);
h34 = zeros(1, NM);
h35 = zeros(1, NM);
s29 = zeros(1, NM);
s30 = zeros(1, NM);
s31 = zeros(1, NM);
s32 = zeros(1, NM);
s33 = zeros(1, NM);
s34 = zeros(1, NM);
s35 = zeros(1, NM);
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
for i=1:NM
    if (i ~= 5)
        mdot29(i) = mdot28(i);
    else
        mdot30(i) = mdot28(i);
    end
end
% Find mass and mole fractions here for Dalton's law
% Mass fractions
for i = 1:NM
    Y29(i) = mdot29(i)/sum(mdot29);
    Y30(i) = mdot30(i)/sum(mdot30);
end
% Mole fractions
X29 = Y2X(Y29, M);
X30 = Y2X(Y30, M);
% Conservation of Momentum
p(29) = p(28);
p(30) = p(28);
% Conservation of Energy
T(29) = T(2);
T(30) = T(2);
% State 29
% Partial pressures have been added as this is a gas mixture.
sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
r29 = RP.ABFLSHdll('TP', T(29), X29(2)*p(29)/1000, zu, int8(0));
h29(2) = (r29.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s29(2) = (r29.s + dsn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
sm = RP.SETFLUIDSdll('METHANE.FLD');
r29 = RP.ABFLSHdll('TP', T(29), X29(9)*p(29)/1000, zu, int8(0));
h29(9) = (r29.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s29(9) = (r29.s + dsn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
sm = RP.SETFLUIDSdll('CO.FLD');
r29 = RP.ABFLSHdll('TP', T(29), X29(6)*p(29)/1000, zu, int8(0));
h29(6) = (r29.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s29(6) = (r29.s + dsn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
% State 30
h30 = h2;
s30 = s2;
% Find the heat transfer
Qdot282930 = 0;
for i=1:NM
    Qdot282930 = Qdot282930 + mdot29(i)*h29(i) + mdot30(i)*h30(i) - mdot28(i)*h28(i);
end
for i=1:NM
    Y29(i) = mdot29(i)/sum(mdot29);
end
X29 = Y2X(Y29, M);

%% Loop here mixing box; There's no mixing box as of 12-03-2021. We removed it.
% Therefore, we have to solve for the compressor directly. However, we do
% not have sufficient information.
%% Depcik addition (12/2/2021) - Let's work backwards
diff = 1;
T(34) = 111;     % This is what we are going to iterate upon
while (diff > 1e-5)
    % Assumptions
    T(34) = T(34)+0.1;
    % What works for T(35): 190.55
    T(32) = TCH4store;     % This is what we are assuming the temperature will be exiting the heat exchanger
    % Find the pressure of liquid methane at the assumed temperature
    % Conservation of momentum
    sm = RP.SETFLUIDSdll('METHANE.FLD');
    r34 = RP.ABFLSHdll('TQ', T(34), 0, zu, int8(0));
    p(34) = r34.P*1000; % [Pa] 
    % Effectively, this would be the partial pressure of methane at State 34
    % Conservation of mass - compressor
    mdot31 = mdot29;
    Y31 = Y29;  % Nothing changes
    X31 = X29;  % Nothing changes
    % Conservation of mass - heat exchanger
    mdot32 = mdot31;
    Y32 = Y31;  % Nothing changes
    X32 = X31;  % Nothing changes
    % Conservation of mass - valve
    mdot33 = mdot32;
    Y33 = Y32;  % Nothing changes
    X33 = X32;  % Nothing changes
    % Find the pressure at State 34
    p(33) = p(34)/X33(9);
    % The temperature at State 34 should be that of State 35 also
    T(33) = T(34);
    % Conservation of Momentum - Vapor-liquid separator
    p(35) = p(33) - p(34);      % Dalton's law of partial pressures
    % Find mdot34*h34 exiting the valve
    sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
    r33 = RP.ABFLSHdll('TP', T(33), X33(2)*p(33)/1000, zu, int8(0));
    h33(2) = (r33.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s33(2) = (r33.s + dsn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('METHANE.FLD');
    r33 = RP.ABFLSHdll('TP', T(33), X33(9)*p(33)/1000, zu, int8(0));
    h33(9) = (r33.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s33(9) = (r33.s + dsn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('CO.FLD');
    r33 = RP.ABFLSHdll('TP', T(33), X33(6)*p(33)/1000, zu, int8(0));
    h33(6) = (r33.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s33(6) = (r33.s + dsn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    m33h33 = 0;
    for i=1:NM
        m33h33 = m33h33 + mdot33(i)*h33(i);
    end
    % Now find the pressure entering the valve
    p(32) = fminsearch(@objectivefcn6, p(33), [], T(32), m33h33, mdot32, M, dhn, NM, X32);
    % Check to make sure there is convergence
    sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
    r32 = RP.ABFLSHdll('TP', T(32), X32(2)*p(32)/1000, zu, int8(0));
    h32(2) = (r32.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s32(2) = (r32.s + dsn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('METHANE.FLD');
    r32 = RP.ABFLSHdll('TP', T(32), X32(9)*p(32)/1000, zu, int8(0));
    h32(9) = (r32.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s32(9) = (r32.s + dsn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('CO.FLD');
    r32 = RP.ABFLSHdll('TP', T(32), X32(6)*p(32)/1000, zu, int8(0));
    h32(6) = (r32.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s32(6) = (r32.s + dsn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    m32h32 = 0;
    for i=1:NM
        m32h32 = m32h32 + mdot32(i)*h32(i);
    end
    diff = abs(m32h32 - m33h33);
    [T(34) diff]
end
T(35) = T(33);
%% SHAH - I STOPPED HERE: THE SUBSEQUENT ITEMS MIGHT NEED TO BE REDONE
%% Taking up the remaining States, i.e., 29 and 31
% State 29 & 31 - compressor
% Conservation of mass
% This is already done inside the loop!
% Conservation of momentum
p(31) = p(32); % [Pa]; We are going in reverse. Since the heat exchanger operates at constant pressure, its inlet pressure is the same as its exit pressure. Also, its inlet is compressor's outlet.
% Conservation of energy
% Need to solve isentropic first and use the compressor efficiency
% Conservation of entropy
m29s29 = 0;
m29h29 = 0;
for i = 1:NM
    m29h29 = m29h29 + mdot29(i)*h29(i); % [W]
    m29s29 = m29s29 + mdot29(i)*s29(i); % [W/K]
end
T31g = T(29)*(p(31)/p(29))^(1.25 - 1); % [K]; Guess temperature
T31s = fminsearch(@objectivefcn2, T31g, [], p(31), m29s29, mdot31, M, dsn, NM, X31);
% Now, we need to find the enthalpy at State 21s
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
h31s = zeros(1, NM);
sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
r31s = RP.ABFLSHdll('TP', T31s, X31(2)*p(31)/1000, zu, int8(0));
h31s(2) = (r31s.h + dhn(2))/M(2)*1000; % [J/kg]
sm = RP.SETFLUIDSdll('CO.FLD');
r31s = RP.ABFLSHdll('TP', T31s, X31(6)*p(31)/1000, zu, int8(0));
h31s(6) = (r31s.h + dhn(6))/M(6)*1000; % [J/kg]
sm = RP.SETFLUIDSdll('METHANE.FLD');
r31s = RP.ABFLSHdll('TP', T31s, X31(9)*p(31)/1000, zu, int8(0));
h31s(9) = (r31s.h + dhn(9))/M(9)*1000; % [J/kg]
m31h31s = 0;
for i = 1:NM
    m31h31s = m31h31s + mdot31(i)*h31s(i); % [W]
end
% Use the isentropic efficiency equation
m31h31 = (m29h29*etac + (m31h31s - m29h29))/etac; % [W]
% Now, we have to find the temperature that achieves this enthalpy
T(31) = fminsearch(@objectivefcn3, T31s, [], p(31), m31h31, mdot31, M, dhn, NM, X31);
% Finally, we have State 32 defined and we can calculate the poewr needed
sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
r31 = RP.ABFLSHdll('TP', T(31), X31(2)*p(31)/1000, zu, int8(0));
h31(2) = (r31.h + dhn(2))/M(2) * 1000; % [J/kg]
s31(2) = (r31.s + dsn(2))/M(2) * 1000; % [J/kg.K]
sm = RP.SETFLUIDSdll('METHANE.FLD');
r31 = RP.ABFLSHdll('TP', T(31), X31(9)*p(31)/1000, zu, int8(0));
h31(9) = (r31.h + dhn(9))/M(9) * 1000; % [J/kg]
s31(9) = (r31.s + dsn(9))/M(9) * 1000; % [J/kg.K]
sm = RP.SETFLUIDSdll('CO.FLD');
r31 = RP.ABFLSHdll('TP', T(31), X31(6)*p(31)/1000, zu, int8(0));
h31(6) = (r31.h + dhn(6))/M(6) * 1000; % [J/kg]
s31(6) = (r31.s + dsn(6))/M(6) * 1000; % [J/kg.K]
Wdot2931 = m29h29 - m31h31; % [W]

% Heat from heat exchanger; States 32 & 33
Qdot3132 = m32h32 - m31h31;

%% Vapor/Liquid separator
% Conservation of mass
mdot34 = zeros(1, NM);
mdot34(9) = mdot33(9);
mdot35 = mdot33 - mdot34;
% Mass fractions
Y34 = zeros(1, NM);
Y35 = zeros(1, NM);
for i = 1:NM
    Y34(i) = mdot34(i)/sum(mdot34);
    Y35(i) = mdot35(i)/sum(mdot35);
end
% Mole fractions
X34 = Y2X(Y34, M);
X35 = Y2X(Y35, M);
% Conservation of momentum
% Already done.
% Conservation of energy
% Already done.
sm = RP.SETFLUIDSdll('METHANE.FLD');
r34 = RP.ABFLSHdll('TP', T(34), p(34), zu, int8(0));
h34(9) = (r34.h + dhn(9))/M(9) * 1000; % [J/kg]
s34(9) = (r34.s + dsn(9))/M(9) * 1000; % [J/kg.K]

sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
r35 = RP.ABFLSHdll('TP', T(35), X35(2) * p(35), zu, int8(0));
h35(2) = (r35.h + dhn(2))/M(2) * 1000; % [J/kg]
s35(2) = (r35.s + dsn(2))/M(2) * 1000; % [J/kg.K]

sm = RP.SETFLUIDSdll('CO.FLD');
r35 = RP.ABFLSHdll('TP', T(35), X35(6) * p(35), zu, int8(0));
h35(6) = (r35.h + dhn(6))/M(6) * 1000; % [J/kg]
s35(6) = (r35.s + dsn(6))/M(6) * 1000; % [J/kg.K]

sigmadot333435 = mdot35*s35' + mdot34*s34' - mdot33*s33';

%% Saving data
mdotall = [mdot1', mdot2', mdot3', mdot4', mdot5', mdot6', mdot7', mdot8', ...
    mdot9', mdot10', mdot11', mdot12', mdot13', mdot14', mdot15', mdot16', ...
    mdot17', mdot18', mdot19', mdot20', mdot21', mdot22', mdot23', mdot24', ...
    mdot25', mdot26', mdot27', mdot28', mdot29', mdot30', mdot31', mdot32', ...
    mdot33', mdot34', mdot35'];

hall = [h1', h2', h3', h4', h5', h6', h7', h8', ...
    h9', h10', h11', h12', h13', h14', h15', h16', ...
    h17', h18', h19', h20', h21', h22', h23', h24', ...
    h25', h26', h27', h28', h29', h30', h31', h32', ...
    h33', h34', h35'];

sall = [s1', s2', s3', s4', s5', s6', s7', s8', ...
    s9', s10', s11', s12', s13', s14', s15', s16', ...
    s17', s18', s19', s20', s21', s22', s23', s24', ...
    s25', s26', s27', s28', s29', s30', s31', s32', ...
    s33', s34', s35'];

%% Add the entropy calculations
% Added by Saud on 12-06

% Which processes need the entropy analysis:
% 
% O2 liquefaction unit [these are already done in the O2 liquefaction
% block]
% 9-10
% 10-11
% 11-12
% 13-14
% 
% Sabatier+
% 7-17-19
% 20-21
% 21-22
% 22-23-24-25
% 24-27
% 15-27-16
% 23-26
% 
% Methane liquefaction unit
% 25-28
% 28-29-30
% 29-31
% 31-32
% 32-33

%% Atmospheric CO2 collection
% sigmadot342 = mdot4*s4' + mdot2*s2' - mdot3*s3';

%% Electrolysis
sigmadot56 = mdot7*s7' + mdot6*s6' - mdot5*s5';

%% Sabatier+
% State 7-26-17: Heat exchanger
sigmadot72617 = mdot17*s17' - mdot26*s26' - mdot7*s7' - Qdot72617/1000;
% State 16-17-18: Sabatier block (rejecting heat to ambient)
sigmadot161718 = mdot18*s18' - mdot17*s17' - mdot16*s16' - Qdot161718/210;
% State 18-19-20: Condenser
sigmadot181920 = mdot20*s20' + mdot19*s19' - mdot18*s18' - Qdot181920/210;
% State 20-21: Compressor after water condensation
sigmadot2021 = mdot21*s21' - mdot20*s20';
% State 21-22: Heat exchanger rejecting heat (perhaps to Martian ambient)
sigmadot2122 = mdot22*s22' - mdot21*s21' - Qdot2122/210;
% State 22-23-24-25: CO2/H2 separator
sigmadot22232425 = mdot25*s25' + mdot24*s24' + mdot23*s23' - mdot22*s22';
% State 24-27: Compressor
sigmadot2427 = mdot27*s27' - mdot24*s24';
% State 15-16-27: Heated valve for carbon dioxide
sigmadot152716 = mdot16*s16' - mdot15*s15' - mdot27*s27' - Qdot152716/1000;
% State 23-26: Valve
sigmadot2326 = mdot26*s26' - mdot23*s23';

%% Methane liquefaction unit
% State 25-28: Compressor
sigmadot2528 = mdot28*s28' - mdot25*s25';
% State 28-29-30: Heat exchanger
sigmadot282930 = mdot30*s30' + mdot29*s29' - mdot28*s28' - Qdot282930/210;
% State 29-31: Compressor
sigmadot2932 = mdot31*s31' - mdot29*s29';
% State 31-32: Heat exchanger
sigmadot3132 = mdot32*s32' - mdot31*s31' - Qdot3132/210;
% State 32-33: Valve
sigmadot3233 = mdot33*s33' - mdot32*s32';