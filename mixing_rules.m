function [a_m, b_m, d_m] = mixing_rules(y_1, y_2, a_ij, b_ij, d_ij, n)

    a_term_1 = (y_1^2)*(a_ij(1,1,n)^(2/3))*(b_ij(1,1,n)^(1/3));
    a_term_2 = 2*y_1*y_2*(a_ij(1,2,n)^(2/3))*(b_ij(1,2,n)^(1/3));
    a_term_3 = (y_2^2)*(a_ij(2,2,n)^(2/3))*(b_ij(2,2,n)^(1/3));
    a_term_4 = sqrt(((y_1^2)*b_ij(1,1,n)) + (2*y_1*y_2*b_ij(1,2,n)) + ((y_2^2)*b_ij(2,2,n)));

    a_m = ((a_term_1 + a_term_2 + a_term_3)^(1.5))/a_term_4;

    b_term_1 = (y_1^2)*b_ij(1,1,n);
    b_term_2 = 2*y_1*y_2*b_ij(1,2,n);
    b_term_3 = (y_2^2)*b_ij(2,2,n);

    b_m = b_term_1 + b_term_2 + b_term_3;

    d_term_1 = (y_1^2)*d_ij(1,1,n);
    d_term_2 = 2*y_1*y_2*d_ij(1,2,n);
    d_term_3 = (y_2^2)*d_ij(2,2,n);

    d_m = d_term_1 + d_term_2 + d_term_3;

end