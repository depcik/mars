function f = objectivefcn5(pCH4, Tliquefy)
%% Initialize REFPROP10
% I don't think partial pressures need to be added here. Instead, this will
% give the partial pressure of methane.
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
sm = RP.SETFLUIDSdll('METHANE.FLD');
rCH4 = RP.ABFLSHdll('TP', Tliquefy, pCH4/1000, zu, int8(0));
f = rCH4.h;