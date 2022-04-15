function f = X2Y(X, M)
% Converts mole fractions (X) to mass fractions (Y).
NM = length(X);
Y = zeros(1, NM);

XiMi = 0;
for i = 1:NM
    XiMi = XiMi + X(i) * M(i);
end

for i = 1:NM
    Y(i) = (X(i)*M(i))/XiMi;
end

f = Y;