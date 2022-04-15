function X = Y2X(Y, M)
% The function Y2X(Y, M, NM) computes the mole fractions 'X' from the mass
% fractions 'Y' using the molecular weights 'M' and the number of species
% specifier 'NM'.

%% Initialize
NM = length(Y);
X = zeros(1, NM);

YiMi = 0;
for i = 1:NM
    YiMi = YiMi + Y(i)/M(i);
end

for i = 1:NM
    X(i) = (Y(i)/M(i))/YiMi;
end