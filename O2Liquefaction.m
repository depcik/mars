function [mdot8, mdot9, mdot10, mdot11, mdot12, mdot13, mdot14, p, T, h8, h9, h10, ...
    h11, h12, h13, h14, s8, s9, s10, s11, s12, s13, s14, q12, Wdot910, Qdot1011, ...
    Qdot1314, X8, X9, X10, X11, X12, X13, X14, Y8, Y9, Y10, Y11, Y12, Y13, Y14, ...
    sigmadot910, sigmadot1011, sigmadot1112, sigmadot1314, sigmadot6149] = O2Liquefaction(NM, M, TO2store, p, T, mdot6, h6, s6, etac, dhn, dsn)
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
% TO2store = temperature [K] at which we want to store the oxygen
% p(8) = pressure [Pa] at which we want to store the oxygen
sm6 = RP.SETFLUIDSdll('OXYGEN.FLD');
T(11) = TO2store;             % Temperature at which we want to cool the oxygen gas down to
p(8) = p(6);
% TMarsNight = 203.15;          % Temperature exiting heat exchangers [K]
TMarsNight = 210; % Setting this to 210 K to match the manuscript values
%% Initialization
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
mdot8 = zeros(1,NM);
mdot9 = zeros(1,NM);
mdot10 = zeros(1,NM);
mdot11 = zeros(1,NM);
mdot12 = zeros(1,NM);
mdot13 = zeros(1,NM);
mdot14 = zeros(1,NM);
Y8 = zeros(1, NM);
Y9 = zeros(1, NM);
Y10 = zeros(1, NM);
Y11 = zeros(1, NM);
Y12 = zeros(1, NM);
Y13 = zeros(1, NM);
Y14 = zeros(1, NM);
h8 = zeros(1,NM);
s8 = zeros(1,NM);
h9 = zeros(1,NM);
s9 = zeros(1,NM);
h10 = zeros(1,NM);
s10 = zeros(1,NM);
h11 = zeros(1,NM);
s11 = zeros(1,NM);
h12 = zeros(1,NM);
s12 = zeros(1,NM);
h13 = zeros(1,NM);
s13 = zeros(1,NM);
h14 = zeros(1,NM);
s14 = zeros(1,NM);
%% Pressure at State 10
% So, now we need to find the pressure that provides the smallest value of
% enthalpy through the valve
x0 = 30*1e6;            % Guess at pressure
p(10) = fminsearch(@objectivefcn1, x0, [], TO2store);        % Pressure to store the oxygen at [Pa]
% Will use p10 = 39 MPa for example calculations
%p(10) = 39*1e6;            % Pa
% First time through we use the mass flow at State 6 to perform our initial
% calcuations. Then, we find the mixing box and re-do all of the calculations
%% Compressor
% First guesses at State 9
T(9) = T(6);   % First guess
mdot9 = mdot6;  % First guess
eps = 1;
% Loop here ****
while (eps > 1e-5)
    % Conservation of Mass
    mdot10 = mdot9;
    % Conservation of Momentum
    p(9) = p(6);
    % Conservation of Energy
    % Steady-state, adiabatic, 0-D flows in and out, no KE, no PE
    r9 = RP.ABFLSHdll('TP', T(9), p(9)/1000, zu, int8(0));
    h9(3) = (r9.h + dhn(3))/M(3)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
    s9(3) = (r9.s + dsn(3))/M(3)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    % Conservation of Entropy
    % Isentropic analysis
    %s10s = s9;
    % Find the isentropic enthalpy
    r10s = RP.ABFLSHdll('PS', p(10)/1000, r9.s, zu, int8(0));     % Pressure [kPa], Entropy [J/molK]
    h10s = (r10s.h + dhn(3))/M(3)*1000;            % [J/kg]
    % Isentropic efficiency
    h10(3) = (h10s-h9(3))/etac + h9(3);     % [J/kg]
    h10m = h10(3)/1000*M(3);                % [J/mol]
    Wdot910 = mdot9(3)*h9(3) - mdot10(3)*h10(3);        % [W]
    h10mR = h10m - dhn(3);  % Taking off the normalization parameter for REFPROP use
    % Find the temperature at State 10
    r10 = RP.ABFLSHdll('PH', p(10)/1000, h10mR, zu, int8(0));        % Pressure [kPa], Enthalpy [J/mol]
    T(10) = r10.T;                      % [K]
    s10(3) = (r10.s + dsn(3))/M(3)*1000;            % [J/kgK]
    % Entropy generation
    sigmadot910 = mdot10(3)*s10(3) - mdot9(3)*s9(3);    % [W/K]
    %% Heat Exchanger
    % Conservation of Mass
    mdot11 = mdot10;
    % Conservation of Momentum
    p(11) = p(10);
    % Conservation of Energy
    % Steady-state, no work, no KE, no PE
    r11 = RP.ABFLSHdll('TP', T(11), p(11)/1000, zu, int8(0));
    h11(3) = (r11.h + dhn(3))/M(3)*1000;            % [J/kg]
    s11(3) = (r11.s + dsn(3))/M(3)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    Qdot1011 = mdot11(3)*h11(3) - mdot10(3)*h10(3);     % [W]
    % Conservation of Entropy
    sigmadot1011 = mdot11(3)*s11(3) - mdot10(3)*s10(3) - Qdot1011/TMarsNight;   % [W/K]
    %% Valve
    % Conservation of Mass
    mdot12 = mdot11;
    % Conservation of Momentum
    p(12) = p(8);
    % Conservation of Energy
    % Steady-state, no heat transfer, no work, no KE, no PE, 1-D flows
    h12(3) = h11(3);
    h12mR = h12(3)/1000*M(3) - dhn(3);
    % Find the quality
    r12 = RP.ABFLSHdll('PH', p(12)/1000, h12mR, zu, int8(0));        % Pressure [kPa], Enthalpy [J/mol]
    q12 = r12.q;    % Vapor quality on a MOLAR basis (moles of vapor/total moles)
    T(12) = r12.T;
    s12(3) = (r12.s + dsn(3))/M(3)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    % Conservation of Entropy
    sigmadot1112 = mdot12(3)*s12(3) - mdot11(3)*s11(3);  % [W/K]
    %% Vapor/Liquid Separator
    p(13) = p(8);
    % Properties at States 8 and 13
    r8 = RP.ABFLSHdll('PQ', p(8)/1000, 0, zu, int8(0));        % Pressure [kPa]
    h8(3) = (r8.h + dhn(3))/M(3)*1000;            % [J/kg]
    s8(3) = (r8.s + dsn(3))/M(3)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    T(8) = r8.T;
    r13 = RP.ABFLSHdll('PQ', p(13)/1000, 1, zu, int8(0));        % Pressure [kPa]
    h13(3) = (r13.h + dhn(3))/M(3)*1000;            % [J/kg]
    s13(3) = (r13.s + dsn(3))/M(3)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    T(13) = r13.T;
    % Find the mass flow rate to the reservoir
    ndot12O2 = mdot12(3)/(M(3)/1000);                       % [kg/s] / ([gm/mol] / [1000 gm/kg]) => [mol/s]
    ndot8O2 = (1-q12)*ndot12O2;
    ndot13O2 = q12*ndot12O2;
    mdot8(3) = ndot8O2*(M(3)/1000);                 % Liquid state
    mdot13(3) = ndot13O2*(M(3)/1000);               % Vapor state
    % Check:
    %chk = mdot13(3) + mdot8(3) - mdot12(3);
    %% Heat Exchanger
    % Conservation of Mass
    mdot14 = mdot13;
    % Conservation of Momentum
    p(14) = p(13);
    T(14) = TMarsNight;             % Temperature at which we want to heat the oxygen gas to
    % Steady-state, no work, no KE, no PE
    r14 = RP.ABFLSHdll('TP', T(14), p(14)/1000, zu, int8(0));
    h14(3) = (r14.h + dhn(3))/M(3)*1000;            % [J/kg]
    s14(3) = (r14.s + dsn(3))/M(3)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
    Qdot1314 = mdot14(3)*h14(3) - mdot13(3)*h13(3);     % [W]
    % Conservation of Entropy
    sigmadot1314 = mdot14(3)*s14(3) - mdot13(3)*s13(3) - Qdot1314/TMarsNight;   % [W/K]
    %% Now, we do our Adiabatic Mixing process (6,14 => 9) to find a new State 9
    % For the system to work in a steady-state scenario, the mass that enters
    % at State 6 must be balanced by the mass that is stored at State 8
    % Therefore, we can do a scaling analysis:
    if (mdot8(3) ~= mdot6(3))
        SC = mdot13(3)/mdot8(3);
        mdot14(3) = SC*mdot6(3);     % [kg/s]
    else
        eps = 0;
    end
    % Conservation of Mass
    mdot9(3) = mdot6(3) + mdot14(3);
    % Conservation of Momentum
    % Already determined
    % Conservation of Energy
    h9(3) = (mdot6(3)*h6(3) + mdot14(3)*h14(3))/mdot9(3);
    h9mR = h9(3)/1000*M(3) - dhn(3);
    r9 = RP.ABFLSHdll('PH', p(12)/1000, h9mR, zu, int8(0));        % Pressure [kPa], Enthalpy [J/mol]
    T(9) = r9.T;
    s9(3) = (r9.s + dsn(3))/M(3)*1000;
    sigmadot6149 = mdot9(3)*s9(3) - mdot14(3)*s14(3) - mdot6(3)*s6(3);  % [W/K]
end

%     Mass and mole fraction calculations
for i = 1:NM
    Y8(i) = mdot8(i)/sum(mdot8);
    Y9(i) = mdot9(i)/sum(mdot9);
    Y10(i) = mdot10(i)/sum(mdot10);
    Y11(i) = mdot11(i)/sum(mdot11);
    Y12(i) = mdot12(i)/sum(mdot12);
    Y13(i) = mdot13(i)/sum(mdot13);
    Y14(i) = mdot14(i)/sum(mdot14);
end
X8 = Y2X(Y8, M);
X9 = Y2X(Y9, M);
X10 = Y2X(Y10, M);
X11 = Y2X(Y11, M);
X12 = Y2X(Y12, M);
X13 = Y2X(Y13, M);
X14 = Y2X(Y14, M);