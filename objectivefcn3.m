function f = objectivefcn3(T21, p21, m21h21, mdot21, M, dhn, NM, X21)
%% Initialize REFPROP10
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
h21 = zeros(1,NM);
%% Find the entropies using partial pressures
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
r21 = RP.ABFLSHdll('TP', T21, X21(2)*p21/1000, zu, int8(0));
h21(2) = (r21.h + dhn(2))/M(2)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
sm = RP.SETFLUIDSdll('CO2.FLD');
r21 = RP.ABFLSHdll('TP', T21, X21(5)*p21/1000, zu, int8(0));
h21(5) = (r21.h + dhn(5))/M(5)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
sm = RP.SETFLUIDSdll('METHANE.FLD');
r21 = RP.ABFLSHdll('TP', T21, X21(9)*p21/1000, zu, int8(0));
h21(9) = (r21.h + dhn(9))/M(9)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
sm = RP.SETFLUIDSdll('CO.FLD');
r21 = RP.ABFLSHdll('TP', T21, X21(6)*p21/1000, zu, int8(0));
h21(6) = (r21.h + dhn(6))/M(6)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
%% Compute m21s21s
m21h21c = 0;
for i=1:NM
    m21h21c = m21h21c + mdot21(i)*h21(i);
end
%% Set the function
f = abs(m21h21 - m21h21c);