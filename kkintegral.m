function delta_n = kkintegral (E,alpha,Eprime)

% reference IEEE Journal of Quantum Electronics, Vol 26, No 1, P113,
% (1990)

loadconstants;
de = 0.001; % increment 1meV
delta_n = zeros(1,length(E));

for j = 1:length(E)
    for i = 1:length(alpha)
        if (E(j)~=Eprime(i))
        delta_n(j) = delta_n(j) + 2*c*hbar/q*alpha(i)/(Eprime(i)^2-E(j)^2)*de*10; 
        % the last factor is to convert alpha in cm-1 to m-1 for calculation
        else
        end
    end
end

end