P = 1:1:300;
T = 308;
for i=1:300
    PR_method_function(T,P(i));
    system('filter_mole_fractions.py');
    fprintf('\033[34mPython process\033[0m at \033[31mP = %d Bar\033[0m and \033[32mT = %d Kelvin\033[0m finished.\n', i, T);
end