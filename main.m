P = 1:1:300;
T = 308;
for i=1:300
    PR_method_function(T,P(i));
    system('filter_mole_fractions.py');
    fprintf('Python process at P = %d Bar and T = %d Kelvin finished.\n', i, T);
end