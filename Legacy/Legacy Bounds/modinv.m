function b = modinv(x,n)
    if gcd(x,n) ~= 1
        error('x has no inverse modulo n');
    end

    [d, a, b] = gcd(x,n);
    b = mod(a,n);
end