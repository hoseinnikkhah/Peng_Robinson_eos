function [true_y_1, initial_y_1] = difference(new_y_1, old_y_1, error)
    diff = new_y_1 - old_y_1;
    diff = abs(diff);  % You had a typo here: 'dif' should be 'diff'
    
    if abs(new_y_1) > abs(old_y_1)
        base = abs(new_y_1);  % Added abs() to ensure positive base
    else
        base = abs(old_y_1);  % Added abs() to ensure positive base
    end
    
    if base == 0
        relative = 0;  % If both values are zero, set relative difference to zero
    else
        relative = (diff/base)*100; 
    end

    if relative <= error
        true_y_1 = new_y_1;
        initial_y_1 = old_y_1;
    else
        true_y_1 = NaN;
        initial_y_1 = NaN;
    end
end