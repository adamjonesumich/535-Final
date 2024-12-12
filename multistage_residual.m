function result = multistage_residual(Isp,e,INDEX,dv_tot,g0,R_guess)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Rs = [0;0;0];
Rs(INDEX) = R_guess;

lambda = @(n) Rs(n)*(1 + e(n)) / (Isp(n) * g0 * ( 2 + Rs(n)*(1-3*e(n)) - Rs(n)^2*(e(n)^2-e(n)) ) );
lambda_to_use = lambda(INDEX);

A = @(n) e(n)^2-e(n);
B = @(n) 1-3*e(n)-(1+e(n))/(lambda_to_use*Isp(n)*g0);

Rlow = @(n) (-B(n) + sqrt(B(n)^2-8*A(n)))/(2*A(n));
Rhigh = @(n) (-B(n) - sqrt(B(n)^2-8*A(n)))/(2*A(n));
R = @(n) max(Rlow(n),Rhigh(n));

Rs(2) = R(2);
Rs(3) = R(3);

result = dv_tot/g0 - Isp(1)*log(Rs(1)) - Isp(2)*log(Rs(2)) - Isp(3)*log(Rs(3));
end