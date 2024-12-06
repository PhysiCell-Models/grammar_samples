function hill_out = hillResponse(signal, x)

base_val = x(1);
min_val = x(2);
max_val = x(3);
decrease_ec50 = x(4);
decrease_pow = x(5);
increase_ec50 = x(6);
increase_pow = x(7);

hill_out = base_val + (max_val-base_val)./(1+(signal./increase_ec50).^(-increase_pow));
hill_out = hill_out + (min_val - hill_out) .*  (signal./decrease_ec50).^decrease_pow ./ (1 + (signal./decrease_ec50).^decrease_pow);

