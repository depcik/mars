function f = objectivefcn1(p10, Tliquefy)
%% Initialize REFPROP10
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
sm6 = RP.SETFLUIDSdll('OXYGEN.FLD');
r8 = RP.ABFLSHdll('TP', Tliquefy, p10/1000, zu, int8(0));
f = r8.h;