function [dhn, dsn] = speciesnorm()
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
%% NIST Standard Information
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
%   https://webbook.nist.gov/
hmfo = [-241.83, 0, 0, 0, -393.52, -110.53, 0, -285.83, -74.6]*1000;                     % NIST Standard Heats of Formation [J/mol]
smfo = [188.84, 130.68, 205.15, 154.84, 213.79, 197.66, 191.61, 69.95, 186.26];           % NIST Standard Entropies of Formation [J/molK]
%% Normalization of Properties
% Let's calculate the normalization parameters so that all working fluids
% have the same reference state. Yeah, you can probably do this in REFPROP,
% but I do not trust it
Tstd = 298.15;          % Standard temperature [K]
pstd = 101325;          % Standard pressure [Pa]
smstd = RP.SETFLUIDSdll('HYDROGEN.FLD');
rstd = RP.ABFLSHdll('TP', Tstd, pstd/1000, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
dhn(2) = hmfo(2) - rstd.h;
dsn(2) = smfo(2) - rstd.s;
smstd = RP.SETFLUIDSdll('OXYGEN.FLD');
rstd = RP.ABFLSHdll('TP', Tstd, pstd/1000, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
dhn(3) = hmfo(3) - rstd.h;
dsn(3) = smfo(3) - rstd.s;
smstd = RP.SETFLUIDSdll('ARGON.FLD');
rstd = RP.ABFLSHdll('TP', Tstd, pstd/1000, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
dhn(4) = hmfo(4) - rstd.h;
dsn(4) = smfo(4) - rstd.s;
smstd = RP.SETFLUIDSdll('CO2.FLD');
rstd = RP.ABFLSHdll('TP', Tstd, pstd/1000, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
dhn(5) = hmfo(5) - rstd.h;
dsn(5) = smfo(5) - rstd.s;
smstd = RP.SETFLUIDSdll('CO.FLD');
rstd = RP.ABFLSHdll('TP', Tstd, pstd/1000, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
dhn(6) = hmfo(6) - rstd.h;
dsn(6) = smfo(6) - rstd.s;
smstd = RP.SETFLUIDSdll('NITROGEN.FLD');
rstd = RP.ABFLSHdll('TP', Tstd, pstd/1000, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
dhn(7) = hmfo(7) - rstd.h;
dsn(7) = smfo(7) - rstd.s;
smstd = RP.SETFLUIDSdll('METHANE.FLD');
rstd = RP.ABFLSHdll('TP', Tstd, pstd/1000, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
dhn(9) = hmfo(9) - rstd.h;
dsn(9) = smfo(9) - rstd.s;
% Water is a bit weird
% Gaseous water
smstd = RP.SETFLUIDSdll('WATER.FLD');
rstd = RP.ABFLSHdll('TQ', Tstd, 1, zu, int8(0));    % Enthalpy returns [J/mol], entropy [J/molK]
dhn(1) = hmfo(1) - rstd.h;
dsn(1) = smfo(1) - rstd.s;
% Liquid water
rstd = RP.ABFLSHdll('TP', Tstd, pstd/1000, zu, int8(0));
dhn(8) = hmfo(8) - rstd.h;
dsn(8) = smfo(8) - rstd.s;