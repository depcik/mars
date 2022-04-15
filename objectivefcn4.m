function f = objectivefcn4(T26, p26, h26, M, dhn)
%% Initialize REFPROP10
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
%% Find the enthalpy gased on the temperature
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
r26 = RP.ABFLSHdll('TP', T26, p26/1000, zu, int8(0));
h26chk = (r26.h + dhn(2))/M(2)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
%% Set the function
f = abs(h26 - h26chk);