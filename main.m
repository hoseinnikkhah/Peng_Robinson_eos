P = 1:1:300;
T = 308;
for i=1:300
    PR_method_function(T,P(i));
    system('filter_mole_fractions.py');
end

