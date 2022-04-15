function f = objectivefcn6(p33, T33, m34h34, mdot33, M, dhn, NM, X33)
%% Initialize REFPROP10
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
h33 = zeros(1,NM);
%% Find the entropies using partial pressures
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
r33 = RP.ABFLSHdll('TP', T33, X33(2)*p33/1000, zu, int8(0));
h33(2) = (r33.h + dhn(2))/M(2)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
sm = RP.SETFLUIDSdll('METHANE.FLD');
r33 = RP.ABFLSHdll('TP', T33, X33(9)*p33/1000, zu, int8(0));
h33(9) = (r33.h + dhn(9))/M(9)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
sm = RP.SETFLUIDSdll('CO.FLD');
r33 = RP.ABFLSHdll('TP', T33, X33(6)*p33/1000, zu, int8(0));
h33(6) = (r33.h + dhn(6))/M(6)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
%% Compute m21s21s
m33h33c = 0;
for i=1:NM
    m33h33c = m33h33c + mdot33(i)*h33(i);
end
%% Set the function
f = abs(m34h34 - m33h33c);