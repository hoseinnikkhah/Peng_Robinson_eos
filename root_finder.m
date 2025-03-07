function [answers] = root_finder()

for p=1:300
    for i=1:mole_length
        A = ABD_308(1,i,p);
        B = ABD_308(2,i,p);
        D = ABD_308(3,i,p);

        C1 = 1;
        C2 = B-1;
        C3 = (D-3*B^2-2*B+A-2*sqrt(A*D));
        C4 = (B^3+B^2-A*B-B*D+2*B*sqrt(A*D));

        coefficients = [C1, C2, C3, C4];
        all_roots = roots(coefficients);

        real_roots = all_roots(abs(imag(all_roots)) < 1e-10);
        real_roots = real(real_roots);

        if length(real_roots) == 1
            Z = real_roots;
        elseif length(real_roots) == 3
            valid_roots = real_roots(real_roots > B);
    
            if isempty(valid_roots)
                Z = max(real_roots);
            elseif length(valid_roots) == 1
                Z = valid_roots;
            else
                Z = max(valid_roots);
            end
        else
            Z = max(real_roots(real_roots > B));
        end

        Z_roots_308(p,i) = Z;

    end
end

