function f = objectivefcn2(T21s, p21, m20s20, mdot21, M, dsn, NM, X21)
%% Initialize REFPROP10
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
s21s = zeros(1,NM);

%% Find the entropies using partial pressures
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
sm = RP.SETFLUIDSdll('HYDROGEN.FLD');
r21s = RP.ABFLSHdll('TP', T21s, X21(2)*p21/1000, zu, int8(0));
s21s(2) = (r21s.s + dsn(2))/M(2)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
sm = RP.SETFLUIDSdll('CO2.FLD');
r21s = RP.ABFLSHdll('TP', T21s, X21(5)*p21/1000, zu, int8(0));
s21s(5) = (r21s.s + dsn(5))/M(5)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
sm = RP.SETFLUIDSdll('METHANE.FLD');
r21s = RP.ABFLSHdll('TP', T21s, X21(9)*p21/1000, zu, int8(0));
s21s(9) = (r21s.s + dsn(9))/M(9)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
sm = RP.SETFLUIDSdll('CO.FLD');
r21s = RP.ABFLSHdll('TP', T21s, X21(6)*p21/1000, zu, int8(0));
s21s(6) = (r21s.s + dsn(6))/M(6)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
%% Compute m21s21s
m21s21s = 0;
for i=1:NM
    m21s21s = m21s21s + mdot21(i)*s21s(i);
end
%% Set the function
f = abs(m20s20 - m21s21s);