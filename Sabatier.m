function [mdot15, mdot16, mdot17, mdot18, mdot19, mdot20, mdot21, mdot22, ...
    mdot23, mdot24, mdot25, mdot26, mdot27, p, T, h15, h16, h17, h18, h19, h20, ...
    h21, h22, h23, h24, h25, h26, h27, s15, s16, s17, s18, s19, s20, s21, s22, s23, ...
    s24, s25, s26, s27, Qdot1516, Qdot717, Qdot161718, Qdot181920, Wdot2021, ...
    Qdot2122, Wdot22232425, Qdot72617, Wdot2427, Qdot152716,...
    Y15, Y16, Y17, Y18, Y19, Y20, Y21, Y22, Y23, Y24, Y25, Y26, Y27,...
    X15, X16, X17, X18, X19, X20, X21, X22, X23, X24, X25, X26, X27] = Sabatier(NM, M, H2Cm, etaCO2, SCH4, pSab, ...
    TSab, etac, etaH2rec, etaCO2rec, p, T, mdot7, h7, h2, s2, dhn, dsn)
%% REFPROP
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
%% Initialization
mdot15 = zeros(1,NM);
mdot16 = zeros(1,NM);
mdot17 = zeros(1,NM);
mdot18 = zeros(1,NM);
mdot19 = zeros(1,NM);
mdot20 = zeros(1,NM);
mdot21 = zeros(1,NM);
mdot22 = zeros(1,NM);
mdot23 = zeros(1,NM);
mdot24 = zeros(1,NM);
mdot25 = zeros(1,NM);
mdot26 = zeros(1,NM);
mdot27 = zeros(1,NM);
Y15 = zeros(1, NM);
Y16 = zeros(1, NM);
Y17 = zeros(1, NM);
Y18 = zeros(1, NM);
Y19 = zeros(1, NM);
Y20 = zeros(1, NM);
Y21 = zeros(1, NM);
Y22 = zeros(1, NM);
Y23 = zeros(1, NM);
Y24 = zeros(1, NM);
Y25 = zeros(1, NM);
Y26 = zeros(1, NM);
Y27 = zeros(1, NM);
X15 = zeros(1, NM);
X16 = zeros(1, NM);
X17 = zeros(1, NM);
X18 = zeros(1, NM);
X19 = zeros(1, NM);
X20 = zeros(1, NM);
X21 = zeros(1, NM);
X22 = zeros(1, NM);
X23 = zeros(1, NM);
X24 = zeros(1, NM);
X25 = zeros(1, NM);
X26 = zeros(1, NM);
X27 = zeros(1, NM);
h15 = zeros(1,NM);
h16 = zeros(1,NM);
h17 = zeros(1,NM);
h18 = zeros(1,NM);
h19 = zeros(1,NM);
h20 = zeros(1,NM);
h21 = zeros(1,NM);
h22 = zeros(1,NM);
h23 = zeros(1,NM);
h24 = zeros(1,NM);
h25 = zeros(1,NM);
h26 = zeros(1,NM);
h27 = zeros(1,NM);
s15 = zeros(1,NM);
s16 = zeros(1,NM);
s17 = zeros(1,NM);
s18 = zeros(1,NM);
s19 = zeros(1,NM);
s20 = zeros(1,NM);
s21 = zeros(1,NM);
s22 = zeros(1,NM);
s23 = zeros(1,NM);
s24 = zeros(1,NM);
s25 = zeros(1,NM);
s26 = zeros(1,NM);
s27 = zeros(1,NM);
%% Initial Calculations
mdoteps = 1e-5;     % Convergence on H2 flow rate to Sabatier process
% Selectivity of CO is assumed to be balance of what was left
SCO = 1-SCH4;           % [-]
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
%% Determine the mass flow rate of CO2 out of the tank
ndotH2in = mdot7(2)/(M(2)/1000);         % [kg/s] / ([gm/mol]*[kg/1000 gm]) = [mol/s]
ndotCO2in = ndotH2in/H2Cm;
mdot15(5) = ndotCO2in*(M(5)/1000);       % [kg/s]
% The properties at State 15 are the same as those at State 2
T(15) = T(2);
p(15) = p(2);
h15 = h2;
s15 = s2;
%% Heat Exchanger for Hydrogen
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
% Conservation of Mass
mdot17 = mdot7;                             % NOTE: this is updated later after the CO2/H2 separator
% Conservation of Momentum
p(17) = p(7);
% Conservation of Energy
T(17) = TSab;
r17 = RP.ABFLSHdll('TP', T(17), p(17)/1000, zu, int8(0));
h17(2) = (r17.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s17(2) = (r17.s + dsn(2))/M(2)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
Qdot717 = mdot17(2)*h17(2) - mdot7(2)*h7(2); % [W]
%% Heated Valve for CO2
sm = RP.SETFLUIDSdll('CO2.FLD');
% Conservation of Mass
mdot16 = mdot15;                                % NOTE: mdot15 is adjusted based on the CO2 recovered from the separator
for i = 1:NM
    Y15(i) = mdot15(i)/sum(mdot15);
    Y16(i) = mdot16(i)/sum(mdot16);
end
X15 = Y2X(Y15, M);
X16 = Y2X(Y16, M);
% Conservation of Momentum
p(16) = p(7);
% Conservation of Energy
T(16) = TSab;
r16 = RP.ABFLSHdll('TP', T(16), p(16)/1000, zu, int8(0));
h16(5) = (r16.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s16(5) = (r16.s + dsn(5))/M(5)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
Qdot1516 = mdot16(5)*h16(5) - mdot15(5)*h15(5); % [W]
%% HERE IS WHERE WE LOOP THE RECOVERED H2/CO2 %%
eps = 100;
'Sabatier Loop'
while (eps > mdoteps)
    mdot17old = mdot17(2);
    %% Sabatier Reaction
    % 4H2 + CO2 => CH4 + 2H2O
    % Conservation of Mass
    % Here we need to use the CO2 efficiency and CH4 selectivity
    ndotCO2out = ndotCO2in - (etaCO2*ndotCO2in);    % [mol/s]
    mdot18(5) = ndotCO2out*(M(5)/1000);             % [kg/s]
    ndotCH4out = SCH4*(ndotCO2in - ndotCO2out);     % [mol/s]
    mdot18(9) = ndotCH4out*(M(9)/1000);             % [kg/s]
    % CO is generated as a side reaction via the Water Gas Shift
    % CO2 + H2 => CO + H2O
    ndotCOout = SCO*(ndotCO2in - ndotCO2out);       % [mol/s]
    mdot18(6) = ndotCOout*(M(6)/1000);              % [kg/s]
    % Now, water can be generated from the Sabatier reaction
    ndotH2Osout = 2*ndotCH4out;
    % and water can be generated during the Water Gas Shift reaction
    ndotH2Owout = ndotCOout;
    % Therefore, the total water generated as a gas
    mdot18(1) = (ndotH2Osout + ndotH2Owout)*(M(1)/1000);    % [kg/s]
    % Finally, we can use the Conservation of Mass to find the H2 not utilized
    mdot18(2) = (mdot17(2) + mdot16(5) - sum(mdot18));          % Matlab native routine used here
    %% Since we have to use Dalton's law of partial pressures, we must calculate the mass and mole fractions right here!
    % Mass fractions
    for i = 1:NM
        Y17(i) = mdot17(i)/sum(mdot17);
        Y18(i) = mdot18(i)/sum(mdot18);
    end
    % Mole fractions
    X17 = Y2X(Y17, M);
    X18 = Y2X(Y18, M);
    % Conservation of Momentum
    p(18) = p(17);
    % Conservation of Energy
    % Let's have the outlet temperature be equal to the inlet temperature and
    % simply reject heat (that could go to the hydrogen heat exchanger and the
    % CO2 heated valve)
    T(18) = T(17);
    %    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
    sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
    r18 = RP.ABFLSHdll('TP', T(18), X18(2)*p(18)/1000, zu, int8(0));
    h18(2) = (r18.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s18(2) = (r18.s + dsn(2))/M(2)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    sm = RP.SETFLUIDSdll('CO2.FLD');
    r18 = RP.ABFLSHdll('TP', T(18), X18(5)*p(18)/1000, zu, int8(0));
    h18(5) = (r18.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s18(5) = (r18.s + dsn(5))/M(5)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    sm = RP.SETFLUIDSdll('METHANE.FLD');
    r18 = RP.ABFLSHdll('TP', T(18), X18(9)*p(18)/1000, zu, int8(0));
    h18(9) = (r18.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s18(9) = (r18.s + dsn(9))/M(9)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    sm = RP.SETFLUIDSdll('CO.FLD');
    r18 = RP.ABFLSHdll('TP', T(18), X18(6)*p(18)/1000, zu, int8(0));
    h18(6) = (r18.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s18(6) = (r18.s + dsn(6))/M(6)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    sm = RP.SETFLUIDSdll('WATER.FLD');
    r18 = RP.ABFLSHdll('TP', T(18), X18(1)*p(18)/1000, zu, int8(0));
    h18(1) = (r18.h + dhn(1))/M(1)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s18(1) = (r18.s + dsn(1))/M(1)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    Qdot161718 = 0;
    for i=1:NM
        Qdot161718 = Qdot161718 + mdot18(i)*h18(i) - mdot16(i)*h16(i) - mdot17(i)*h17(i);
    end
%     % Mass fractions exiting - THIS IS NOT NEEDED HERE SINCE WE HAVE
%     CALCULATED THIS EARLIER! 
%     for i=1:NM
%         Y18(i) = mdot18(i)/sum(mdot18);
%     end
    %% Water Recovery
    % It makes sense that we want to remove the water prior to separating out
    % the CO2 and H2 since we want dry streams of those. Moreover, we need to
    % cool down the mixture for the CO2/H2 separator
    % Conservation of Mass
    %    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
    for i=1:NM
        if (i ~= 1)
            mdot20(i) = mdot18(i);
        else
            mdot19(8) = mdot18(1);      % We are condensing the water
        end
    end
    %% Again, Dalton's law of partial pressures requires that we compute mass and mole fractions right here!
    % Mass fractions
    for i = 1:NM
        Y19(i) = mdot19(i)/sum(mdot19);
        Y20(i) = mdot20(i)/sum(mdot20);
    end
    % Mole fractions
    X19 = Y2X(Y19, M);
    X20 = Y2X(Y20, M);
    % Conservation of Momentum
    p(19) = p(1);
    p(20) = p(1);
    % Conservation of Energy
    T(19) = T(1);
    T(20) = T(1);
    %    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
    sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
    r20 = RP.ABFLSHdll('TP', T(20), X20(2)*p(20)/1000, zu, int8(0));
    h20(2) = (r20.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s20(2) = (r20.s + dsn(2))/M(2)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    sm = RP.SETFLUIDSdll('CO2.FLD');
    r20 = RP.ABFLSHdll('TP', T(20), X20(5)*p(20)/1000, zu, int8(0));
    h20(5) = (r20.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s20(5) = (r20.s + dsn(5))/M(5)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    sm = RP.SETFLUIDSdll('METHANE.FLD');
    r20 = RP.ABFLSHdll('TP', T(20), X20(9)*p(20)/1000, zu, int8(0));
    h20(9) = (r20.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s20(9) = (r20.s + dsn(9))/M(9)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    sm = RP.SETFLUIDSdll('CO.FLD');
    r20 = RP.ABFLSHdll('TP', T(20), X20(6)*p(20)/1000, zu, int8(0));
    h20(6) = (r20.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s20(6) = (r20.s + dsn(6))/M(6)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    sm = RP.SETFLUIDSdll('WATER.FLD');
    r19 = RP.ABFLSHdll('TP', T(19), p(19)/1000, zu, int8(0)); % WE DON'T NEED TO USE PARTIAL PRESSURES HERE, RIGHT?
    h19(8) = (r19.h + dhn(8))/M(8)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s19(8) = (r19.s + dsn(8))/M(8)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    Qdot181920 = 0;
    for i=1:NM
        Qdot181920 = Qdot181920 + mdot20(i)*h20(i) + mdot19(i)*h19(i) - mdot18(i)*h18(i);
    end
    %% Compressor before CO2/H2 Separator
    % Conservation of Mass
    mdot21 = mdot20;
    % Dalton's law of partial pressures
    % Mass fractions 
    for i = 1:NM
        Y21(i) = mdot21(i)/sum(mdot21);
    end
    % Mole fractions
    X21 = Y2X(Y21, M);
    % Conservation of Momentum
    p(21) = p(22);
    % Conservation of Energy
    % Need to solve isentropic first and use the compressor efficiency
    % Conservation of Entropy
    m20s20 = 0;
    m20h20 = 0;
    for i=1:NM
        m20s20 = m20s20 + mdot20(i)*s20(i);     % [W/K]
        m20h20 = m20h20 + mdot20(i)*h20(i);     % [W]
    end
    T21g = T(20)*((p(22)/p(20))^(1.25-1));  % Guess of isentropic temperature: State 21s [K]
    T21s = fminsearch(@objectivefcn2, T21g, [], p(21), m20s20, mdot21, M, dsn, NM, X21);
    % Now, we need to find the enthalpy at State 21s
    %    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
    % Use mole fractions to add the Dalton's law functionality
    h21s = zeros(1,NM);
    sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
    r21s = RP.ABFLSHdll('TP', T21s, X21(2)*p(21)/1000, zu, int8(0));
    h21s(2) = (r21s.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('CO2.FLD');
    r21s = RP.ABFLSHdll('TP', T21s, X21(5)*p(21)/1000, zu, int8(0));
    h21s(5) = (r21s.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('METHANE.FLD');
    r21s = RP.ABFLSHdll('TP', T21s, X21(9)*p(21)/1000, zu, int8(0));
    h21s(9) = (r21s.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('CO.FLD');
    r21s = RP.ABFLSHdll('TP', T21s, X21(6)*p(21)/1000, zu, int8(0));
    h21s(6) = (r21s.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    m21h21s = 0;
    for i=1:NM
        m21h21s = m21h21s + mdot21(i)*h21s(i);     % [W]
    end
    % Use the isentropic efficiency equation
    m21h21 = (m20h20*etac + (m21h21s-m20h20))/etac;   % [W]
    % Now, we have to find a temperature that achieves this enthalpy
    T(21) = fminsearch(@objectivefcn3, T21s, [], p(21), m21h21, mdot21, M, dhn, NM, X21);
    % Finally, we have State 21 defined and we can calculate the power needed
    % Partial pressure functionality has been added to the lines 296 to
    % 312.
    sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
    r21 = RP.ABFLSHdll('TP', T(21), X21(2)*p(21)/1000, zu, int8(0));
    h21(2) = (r21.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s21(2) = (r21.s + dsn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('CO2.FLD');
    r21 = RP.ABFLSHdll('TP', T(21), X21(5)*p(21)/1000, zu, int8(0));
    h21(5) = (r21.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s21(5) = (r21.s + dsn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('METHANE.FLD');
    r21 = RP.ABFLSHdll('TP', T(21), X21(9)*p(21)/1000, zu, int8(0));
    h21(9) = (r21.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s21(9) = (r21.s + dsn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('CO.FLD');
    r21 = RP.ABFLSHdll('TP', T(21), X21(6)*p(21)/1000, zu, int8(0));
    h21(6) = (r21.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s21(6) = (r21.s + dsn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    Wdot2021 = m20h20 - m21h21;
    %% Heat Exchanger before H2/O2 Separator
    % We want to get the input similar to Sircar & Golden given their efficiencies
    % Conservation of Mass
    mdot22 = mdot21;
    % Mass fractions
    for i = 1:NM
        Y22(i) = mdot22(i)/sum(mdot22);
    end
    % Mole fractions
    X22 = Y2X(Y22, M);
    % Conservation of Momentum
    % p(22) defined prior
    % Conservation of Energy
    % T(22) defined prior
    % Partial pressure functionality has been added here.
    sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
    r22 = RP.ABFLSHdll('TP', T(22), X22(2)*p(22)/1000, zu, int8(0));
    h22(2) = (r22.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s22(2) = (r22.s + dsn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('CO2.FLD');
    r22 = RP.ABFLSHdll('TP', T(22), X22(5)*p(22)/1000, zu, int8(0));
    h22(5) = (r22.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s22(5) = (r22.s + dsn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('METHANE.FLD');
    r22 = RP.ABFLSHdll('TP', T(22), X22(9)*p(22)/1000, zu, int8(0));
    h22(9) = (r22.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s22(9) = (r22.s + dsn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('CO.FLD');
    r22 = RP.ABFLSHdll('TP', T(22), X22(6)*p(22)/1000, zu, int8(0));
    h22(6) = (r22.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s22(6) = (r22.s + dsn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    m22h22 = 0;
    for i=1:NM
        m22h22 = m22h22 + mdot22(i)*h22(i);     % [W]
    end
    % Find the heat transfer
    Qdot2122 = m22h22 - m21h21; % [W]
    %% Separator
    %    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
    % Conservation of Mass
    mdot23(2) = etaH2rec*mdot22(2);
    mdot24(5) = etaCO2rec*mdot22(5);
    mdot25(2) = mdot22(2) - mdot23(2);
    mdot25(5) = mdot22(5) - mdot24(5);
    mdot25(9) = mdot22(9);
    mdot25(6) = mdot22(6);
    % Mass fractions
    for i = 1:NM
        Y23(i) = mdot23(i)/sum(mdot23);
        Y24(i) = mdot24(i)/sum(mdot24);
        Y25(i) = mdot25(i)/sum(mdot25);
    end
    % Mole fractions
    X23 = Y2X(Y23, M);
    X24 = Y2X(Y24, M);
    X25 = Y2X(Y25, M);
    % Conservation of Momentum
    p(23) = p(22);
    p(24) = 101325;     % Pa
    p(25) = 101325;     % Pa
    % Conservation of Energy
    % All of these temperatures are assumed
    T(23) = T(22);
    T(24) = T(22);
    T(25) = T(22);
    sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
    r23 = RP.ABFLSHdll('TP', T(23), p(23)/1000, zu, int8(0)); % PARTIAL PRESSURE NOT USED HERE BECAUSE SINGLE STREAM!
    h23(2) = (r23.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s23(2) = (r23.s + dsn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    r25 = RP.ABFLSHdll('TP', T(25), X25(2)*p(25)/1000, zu, int8(0)); % PARTIAL PRESSURES USED HERE!
    h25(2) = (r25.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s25(2) = (r25.s + dsn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('CO2.FLD');
    r24 = RP.ABFLSHdll('TP', T(24), p(24)/1000, zu, int8(0)); % PARTIAL PRESSURE NOT USED HERE BECAUSE SINGLE STREAM!
    h24(5) = (r24.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s24(5) = (r24.s + dsn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    r25 = RP.ABFLSHdll('TP', T(25), X25(5)*p(25)/1000, zu, int8(0)); % PARTIAL PRESSURES USED HERE!
    h25(5) = (r25.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s25(5) = (r25.s + dsn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('METHANE.FLD');
    r25 = RP.ABFLSHdll('TP', T(25), X25(9)*p(25)/1000, zu, int8(0)); % PARTIAL PRESSURES USED HERE!
    h25(9) = (r25.h + dhn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s25(9) = (r25.s + dsn(9))/M(9)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    sm = RP.SETFLUIDSdll('CO.FLD');
    r25 = RP.ABFLSHdll('TP', T(25), X25(6)*p(25)/1000, zu, int8(0)); % PARTIAL PRESSURES USED HERE!
    h25(6) = (r25.h + dhn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s25(6) = (r25.s + dsn(6))/M(6)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    % Compute the work it takes to get everything to these states
    Wdot22232425 = 0;
    for i=1:NM
        Wdot22232425 = Wdot22232425 + mdot22(i)*h22(i) - mdot23(i)*h23(i) - mdot24(i)*h24(i) - mdot25(i)*h25(i);
    end
    %% Hydrogen valve
    %    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
    % Conservation of Mass
    mdot26 = mdot23;
    % Mass fractions
    for i = 1:NM
        Y26(i) = mdot26(i)/sum(mdot26);
    end
    % Mole fractions
    X26 = Y2X(Y26, M);
    % Conservation of Momentum
    p(26) = p(7);
    % Conservation of Energy
    % Adiabatic valve
    h26(2) = h23(2);
    T26g = T(7);
    % Find the temperature that equates to this enthalpy given the pressure
    T(26) = fminsearch(@objectivefcn4, T26g, [], p(26), h26(2), M, dhn);
    % Check
    sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
    r26 = RP.ABFLSHdll('TP', T(26), p(26)/1000, zu, int8(0));
    h26(2) = (r26.h + dhn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s26(2) = (r26.s + dsn(2))/M(2)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    %% Hydrogen mixing (7, 26 => 17)
    % Conservation of Mass
    mdot17(2) = mdot26(2) + mdot7(2);
    % Mass fractions
    for i = 1:NM
        Y17(i) = mdot17(i)/sum(mdot17);
    end
    % Mole fractions
    X17 = Y2X(Y17, M);
    % Conservation of Momentum
    % All are at the same pressure
    % Conservation of Energy
    % 0 = Qdot72617 + mdot7*h7 + mdot26*h26 - mdot17*h17
    % State 17 is set in terms of temperature and pressure, so that does not
    % change. Hence, the value of enthalpy that we already calculated is good
    % to go
    Qdot72617 = mdot17(2)*h17(2) - mdot7(2)*h7(2) - mdot26(2)*h26(2);       % [W]
    % LOOP NEEDED - this new value of mdot17 will cause everything to change
    % AND we need to adjust the mass flow rate of the CO2 entering
    %% CO2 Compression (24 => 27)
    %    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
    % Conservation of Mass
    mdot27 = mdot24;
    % Mass fraction
    for i = 1:NM
        Y27(i) = mdot27(i)/sum(mdot27);
    end
    % Mole fraction
    X27 = Y2X(Y27, M);
    % Conservation of Momentum
    p(27) = p(15);
    % Conservation of Energy
    % Need to solve isentropic first and use the compressor efficiency
    s27s = s24(5);      % Isentropic
    s27ms = s27s/1000*M(5) - dsn(5);        % [J/molK]
    % We can use this and the pressure to find the isentropic enthalpy
    sm = RP.SETFLUIDSdll('CO2.FLD');
    r27s = RP.ABFLSHdll('PS', p(27)/1000, s27ms, zu, int8(0));         % Pressure in [kPa], Entropy in [J/molK]
    h27s = (r27s.h + dhn(5))/M(5)*1000;      % [J/kg]
    % Use the isentropic efficiency equation to find the actual enthalpy
    h27(5) = (h24(5)*etac + (h27s-h24(5)))/etac;    % Enthalpy [J/kg]
    h27m = h27(5)/1000*M(5) - dhn(5);         % [J/mol]
    r27 = RP.ABFLSHdll('PH', p(27)/1000, h27m, zu, int8(0));         % Pressure in [kPa], Enthalpy in [J/mol]
    h27(5) = (r27.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s27(5) = (r27.s + dsn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    T(27) = r27.T;
    % Find the work required
    Wdot2427 = mdot24(5)*h24(5) - mdot27(5)*h27(5);     % [W]
    %% Need to find new mdot15 while taking into account mdot27
    ndotH2in = mdot17(2)/(M(2)/1000);           % Molar flow rate of hydrogen taking into account H2 recovered
    ndotCO2in = ndotH2in/H2Cm;                  % The molar flow rate of CO2 required to balance the molar H2 flow rate [mol/s]
    mdot16(5) = ndotCO2in*(M(5)/1000);           % The total mass flow rate into the Sabatier process
    %% CO2 Mixing and Heated Valve redo (15, 27 => 16)
    % Conservation of Mass
    mdot15(5) = mdot16(5) - mdot27(5);          % [kg/s] - the amount that comes from the CO2 tank
    % The following mass and mole fractions were first calculated outside
    % the loop; however, they need to be updated here so that the While
    % loop can function properly. Also, having these here will spit out the
    % updated values for Y and X at States 15 and 16.
    for i = 1:NM
        Y15(i) = mdot15(i)/sum(mdot15);
        Y16(i) = mdot16(i)/sum(mdot16);
    end
    X15 = Y2X(Y15, M);
    X16 = Y2X(Y16, M);

    % Conservation of Momentum
    % Already done
    % Conservation of Energy
    % 0 = Qdot152716 + mdot15*h15 + mdot27*h27 - mdot16*h16
    % State 16 is already defined, so we do not need to redo it
    Qdot152716 = mdot16(5)*h16(5) - mdot15(5)*h15(5) - mdot27(5)*h27(5);    % [W]
    %% LOOP BACK HERE
    eps = abs(mdot17(2) - mdot17old)
end