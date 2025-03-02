y_1 = 2.00807*10^-6;
y_2 = 1 - y_1;
run('PR_Mansoori.m');
% only at temp 308 and pressure 120
T_fixed = 308;
P_fixed = 120;

upper = (y_1^2 * a_ij_3D(1,1,1)^(2/3) * b_ij_3D(1,1,1)^(1/3)) + ...
        (2*y_1*y_2 * a_ij_3D(1,2,1)^(2/3) * b_ij_3D(1,2,1)^(1/3)) + ...
        (y_2^2 * a_ij_3D(2,2,1)^(2/3) * b_ij_3D(2,2,1)^(1/3));
        lower = (y_1^2 * b_ij_3D(1,1,1)) + (2*y_1*y_2*b_ij_3D(1,2,1)) + (y_2^2 * b_ij_3D(2,2,1));

        a_m = (upper^1.5) / sqrt(lower);
        b_m = (y_1^2 * b_ij_3D(1,1,1)) + (2*y_1*y_2*b_ij_3D(1,2,1)) + (y_2^2 * b_ij_3D(2,2,1));
        d_m = (y_1^2 * d_ij_3D(1,1,1)) + (2*y_1*y_2*d_ij_3D(1,2,1)) + (y_2^2 * d_ij_3D(2,2,1));

[A, B, D] = correlations_cap(T_fixed, P_fixed, a_m, b_m, d_m);


for i=1:24
    A_val = A;
    B_val = B;
    D_val = D;
    
    c3 = 1;
    c2 = B_val - 1;
    c1 = D_val - 3*B_val^2 - 2*B_val + A_val - 2*sqrt(A_val*D_val);
    c0 = B_val^3 + B_val^2 - A_val*B_val - B_val*D_val + 2*B_val*sqrt(A_val*D_val);
    
    coeff = [c3, c2, c1, c0];
    root_values = roots(coeff);
    
    % Select only real roots
    real_roots = root_values(imag(root_values) == 0);
    
    % For gas phase calculations, typically the largest real root is used
    if ~isempty(real_roots)
        Z(1,i) = max(real_roots);
    else
        Z(1,i) = NaN; % Handle case where no real roots exist
    end
end

v = (real_roots*R*T_fixed)/P_fixed;

term_1 = y_1*(b_ij_3D(1,1,1) + b_ij_3D(2,1,1)) + y_2*(b_ij_3D(1,2,1) + b_ij_3D(2,2,1));
term_2 = y_1*(a_ij_3D(1,1,1) + a_ij_3D(2,1,1)) + y_2*(a_ij_3D(1,2,1) + a_ij_3D(2,2,1));
term_3 = y_1*(d_ij_3D(1,1,1) + d_ij_3D(2,1,1)) + y_2*(d_ij_3D(1,2,1) + d_ij_3D(2,2,1));

tau = (a_m + R*T_fixed*d_m) - (2*sqrt(a_m*d_m*R*T_fixed)*(1/2 - (term_1/b_m))) + term_2*(1 - sqrt((R*T_fixed*d_m)/a_m)) + term_3*(R*T_fixed - sqrt((R*T_fixed*a_m)/d_m));

coeff1 = 2*(y_1*b_ij_3D(1,1,1) + y_2*b_ij_3D(1,2,1));
coeff2 = (real_roots + (1 + sqrt(2))*B)/(real_roots + (1 - sqrt(2))*B);
        
phi = exp((real_roots-1)*((coeff1/b_m) - 1) - log(real_roots - B) - (tau/(sqrt(2)*R*T_fixed*b_m))*log(coeff2));