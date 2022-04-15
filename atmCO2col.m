function [h2, s2, mdot3, h3, s3, mdot4, h4, s4, p4, T4, X2, X3, X4, Y2, Y3, Y4] = atmCO2col(NM, M, Xatm, mdot2, T, p, dhn, dsn, eta2)
%% Initialize REFPROP10
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Program Files (x86)\REFPROP');
MASSI = RP.MASS_BASE_SI;        % Calculations in Mass basis
iMass = int8(0);                % 0: molar fractions; 1: mass fractions
iFlag = int8(0);                % 0: don't call SATSPLN; 1: call SATSPLN
ierr = int8(0);
zu = {1.0};                     % Composition
sm = RP.SETFLUIDSdll('CO2.FLD');
%% Set up system
h2 = zeros(1,NM);
s2 = zeros(1,NM);
mdot3 = zeros(1,NM);
mdot4 = zeros(1,NM);
h3 = zeros(1,NM);
s3 = zeros(1,NM);
h4 = zeros(1,NM);
s4 = zeros(1,NM);
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
% Check the CO2 in the tank to make sure it is a liquid
r2 = RP.ABFLSHdll('TP', T(2), p(2)/1000, zu, int8(0));
h2(5) = (r2.h + dhn(5))/M(5)*1000;            % [J/mol] / [gm/mol] * [1000 gm/kg] => [J/kg]
s2(5) = (r2.s + dsn(5))/M(5)*1000;            % [J/molK] / [gm/mol] * [1000 gm/kg] => [J/kgK]
if (r2.q == -998)
    'CO2 Tank - Subcooled Liquid'          % It should be a subcooled liquid for storage
end
% Find the exiting gas from the capture process
Mbtm = 0;
for i=1:length(Xatm)
    Mbtm = Mbtm + Xatm(i)*M(i);
end
% Atmospheric mass fractions
Yatm = zeros(1, NM);
for i=1:length(Xatm)
    Yatm(i) = Xatm(i)*M(i)/Mbtm;
end
%    H2O(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7)
% OK, if etac is the collection efficiency and mdot2 is the mass flow rate
% of CO2 that is heading to the tank. mdot3 is the mass flow rate into
% the tank, and mdot4 is that exiting back to the atmosphere
% mdot3(5) = mdot2(5) + mdot4(5);
% eta2 = mdot2(5)/mdot3(5);     % CO2 capture efficiency
% eta2 = mdot2(5)/(mdot2(5) + mdot4(5)
mdot4(5) = (mdot2(5) - eta2*mdot2(5))/eta2; % Mass flow rate of CO2 back to atmosphere [kg/s]
mdot3(5) = mdot2(5) + mdot4(5);     % Mass flow rate of CO2 in [kg/s]
mdot3ttl = mdot3(5)/Yatm(5);        % Total mass flow rate in [kg/s]
mdot4ttl = mdot3ttl - mdot2(5);     % Total mass flow rate out [kg/s]
% Now we can find the rest of the inlet flow rates
Y4 = zeros(1, NM);
X4btm = 0;
for i=1:length(Xatm)
    mdot3(i) = mdot3ttl*Yatm(i);
    if (i ~= 5)
        mdot4(i) = mdot3(i);
    end
    Y4(i) = mdot4(i)/mdot4ttl;
    X4btm = X4btm + Y4(i)/M(i);
end
% Find the exit mole fractions
X4 = zeros(1, NM);
for i=1:length(Xatm)
    X4(i) = (Y4(i)/M(i))/X4btm;
end

% Find mole fractions at 2 and 3 also
Y2 = zeros(1, NM);
Y3 = zeros(1, NM);
for i = 1:NM
    Y2(i) = mdot2(i)/sum(mdot2);
    Y3(i) = mdot3(i)/sum(mdot3);
end

Y2M2 = 0;
Y3M3 = 0;
for i = 1:NM
    Y2M2 = Y2M2 + Y2(i)/M(i);
    Y3M3 = Y3M3 + Y3(i)/M(i);
end

X2 = zeros(1, NM);
X3 = zeros(1, NM);
for i = 1:NM
    X2(i) = (Y2(i)/M(i))/Y2M2;
    X3(i) = (Y3(i)/M(i))/Y3M3;
end

% The pressure and temperature is assumed constant for the process
p4 = p(3);
T4 = T(3);

%% Let's find the properties at State 3 and State 4
% I have tried numerous ways to use REFPROP and generate mixtures, but it
% seems buggy and I haven't had much sustained success. Moreover, REFPROP
% reference states are not well described and tend to be all over the place
% If you want mixture properties, I suggest you do it the long way. This
% also allows you to check the properties of each component separately
%    H2Og(1), H2(2), O2(3), Ar(4), CO2(5), CO(6), N2(7), H2Ol(8), CH4(9)
% NOTE: I only did a few of the species; one could set this up to calculate the values for all species

%% PARTIAL PRESSURES ARE USED HERE TO FIND THE PROPERTIES OF THESE GASES AT THE SPECIFIED STATES
sm = RP.SETFLUIDSdll('CO2.FLD');
r3 = RP.ABFLSHdll('TP', T(3), X3(5)*p(3)/1000, zu, int8(0));
h3(5) = (r3.h + dhn(5))/M(5)*1000;            % [J/kg]
s3(5) = (r3.s + dsn(5))/M(5)*1000;            % [J/kgK]
r4 = RP.ABFLSHdll('TP', T4, X4(5)*p4/1000, zu, int8(0));
h4(5) = (r4.h + dhn(5))/M(5)*1000;            % [J/kg]
s4(5) = (r4.s + dsn(5))/M(5)*1000;            % [J/kgK]
sm = RP.SETFLUIDSdll('NITROGEN.FLD');
r3 = RP.ABFLSHdll('TP', T(3), X3(7)*p(3)/1000, zu, int8(0));
h3(7) = (r3.h + dhn(7))/M(7)*1000;            % [J/kg]
s3(7) = (r3.s + dsn(7))/M(7)*1000;            % [J/kgK]
r4 = RP.ABFLSHdll('TP', T4, X4(7)*p4/1000, zu, int8(0));
h4(7) = (r4.h + dhn(7))/M(7)*1000;            % [J/kg]
s4(7) = (r4.s + dsn(7))/M(7)*1000;            % [J/kgK]
sm = RP.SETFLUIDSdll('ARGON.FLD');
r3 = RP.ABFLSHdll('TP', T(3), X3(4)*p(3)/1000, zu, int8(0));
h3(4) = (r3.h + dhn(4))/M(4)*1000;            % [J/kg]
s3(4) = (r3.s + dsn(4))/M(4)*1000;            % [J/kgK]
r4 = RP.ABFLSHdll('TP', T4, X4(4)*p4/1000, zu, int8(0));
h4(4) = (r4.h + dhn(4))/M(4)*1000;            % [J/kg]
s4(4) = (r4.s + dsn(4))/M(4)*1000;            % [J/kgK]
sm = RP.SETFLUIDSdll('OXYGEN.FLD');
r3 = RP.ABFLSHdll('TP', T(3), X3(3)*p(3)/1000, zu, int8(0));
h3(3) = (r3.h + dhn(3))/M(3)*1000;            % [J/kg]
s3(3) = (r3.s + dsn(3))/M(3)*1000;            % [J/kgK]
r4 = RP.ABFLSHdll('TP', T4, X4(3)*p4/1000, zu, int8(0));
h4(3) = (r4.h + dhn(3))/M(3)*1000;            % [J/kg]
s4(3) = (r4.s + dsn(3))/M(3)*1000;            % [J/kgK]
sm = RP.SETFLUIDSdll('CO.FLD');
r3 = RP.ABFLSHdll('TP', T(3), X3(6)*p(3)/1000, zu, int8(0));
h3(6) = (r3.h + dhn(6))/M(6)*1000;            % [J/kg]
s3(6) = (r3.s + dsn(6))/M(6)*1000;            % [J/kgK]
r4 = RP.ABFLSHdll('TP', T4, X4(6)*p4/1000, zu, int8(0));
h4(6) = (r4.h + dhn(6))/M(6)*1000;            % [J/kg]
s4(6) = (r4.s + dsn(6))/M(6)*1000;            % [J/kgK]


%-------------------------------------------------------------------------
% ORIGINAL CODE BELOW THIS LINE
% sm = RP.SETFLUIDSdll('CO2.FLD');
% r3 = RP.ABFLSHdll('TP', T(3), p(3)/1000, zu, int8(0));
% h3(5) = (r3.h + dhn(5))/M(5)*1000;            % [J/kg]
% s3(5) = (r3.s + dsn(5))/M(5)*1000;            % [J/kgK]
% r4 = RP.ABFLSHdll('TP', T4, p4/1000, zu, int8(0));
% h4(5) = (r4.h + dhn(5))/M(5)*1000;            % [J/kg]
% s4(5) = (r4.s + dsn(5))/M(5)*1000;            % [J/kgK]
% sm = RP.SETFLUIDSdll('NITROGEN.FLD');
% r3 = RP.ABFLSHdll('TP', T(3), p(3)/1000, zu, int8(0));
% h3(7) = (r3.h + dhn(7))/M(7)*1000;            % [J/kg]
% s3(7) = (r3.s + dsn(7))/M(7)*1000;            % [J/kgK]
% r4 = RP.ABFLSHdll('TP', T4, p4/1000, zu, int8(0));
% h4(7) = (r4.h + dhn(7))/M(7)*1000;            % [J/kg]
% s4(7) = (r4.s + dsn(7))/M(7)*1000;            % [J/kgK]
% sm = RP.SETFLUIDSdll('ARGON.FLD');
% r3 = RP.ABFLSHdll('TP', T(3), p(3)/1000, zu, int8(0));
% h3(4) = (r3.h + dhn(4))/M(4)*1000;            % [J/kg]
% s3(4) = (r3.s + dsn(4))/M(4)*1000;            % [J/kgK]
% r4 = RP.ABFLSHdll('TP', T4, p4/1000, zu, int8(0));
% h4(4) = (r4.h + dhn(4))/M(4)*1000;            % [J/kg]
% s4(4) = (r4.s + dsn(4))/M(4)*1000;            % [J/kgK]
% sm = RP.SETFLUIDSdll('OXYGEN.FLD');
% r3 = RP.ABFLSHdll('TP', T(3), p(3)/1000, zu, int8(0));
% h3(3) = (r3.h + dhn(3))/M(3)*1000;            % [J/kg]
% s3(3) = (r3.s + dsn(3))/M(3)*1000;            % [J/kgK]
% r4 = RP.ABFLSHdll('TP', T4, p4/1000, zu, int8(0));
% h4(3) = (r4.h + dhn(3))/M(3)*1000;            % [J/kg]
% s4(3) = (r4.s + dsn(3))/M(3)*1000;            % [J/kgK]
% sm = RP.SETFLUIDSdll('CO.FLD');
% r3 = RP.ABFLSHdll('TP', T(3), p(3)/1000, zu, int8(0));
% h3(6) = (r3.h + dhn(6))/M(6)*1000;            % [J/kg]
% s3(6) = (r3.s + dsn(6))/M(6)*1000;            % [J/kgK]
% r4 = RP.ABFLSHdll('TP', T4, p4/1000, zu, int8(0));
% h4(6) = (r4.h + dhn(6))/M(6)*1000;            % [J/kg]
% s4(6) = (r4.s + dsn(6))/M(6)*1000;            % [J/kgK]