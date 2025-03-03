P = 1:1:300;
T = 308;
for i=1:300
    PR_method_function(T,P(i));
    system('filter_mole_fractions.py');
    fprintf('<a href="matlab:" style="color:blue">Python process</a> at <a href="matlab:" style="color:red">P = %d Bar</a> and <a href="matlab:" style="color:darkgreen">T = %d Kelvin</a> finished.\n', i, T);
end