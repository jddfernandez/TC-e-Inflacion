% Definimos una funcion MCO generica que permite calcular los coeficientes 
% de manera matricial junto a los residuos del modelo

function [beta_gorro, e_gorro,s2,k] = OLS(Y,X)

% Por formula, beta = (X'X)^(-1) (X'Y)
beta_gorro = (X' * X)\(X' * Y);

% Por formula, e = Y - X * beta_gorro
e_gorro = Y - (X * beta_gorro);

k=size(X,2);
n=length(Y);
s2=(e_gorro'*e_gorro)/(n-k);

end  