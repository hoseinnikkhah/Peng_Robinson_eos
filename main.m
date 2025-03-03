P = 1:1:300;
T = 308;
for i=1:300
    PR_method_function(T,P(i));
    system('filter_mole_fractions.py');
    disp(['<strong><font color="blue">Python process</font></strong> at <strong><font color="red">P = ' num2str(i) ' Bar</font></strong> and <strong><font color="darkgreen">T = ' num2str(T) ' Kelvin</font></strong> finished.']);
end