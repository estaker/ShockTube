function A = AFlux(v)
A = zeros(size(v,3));
A(1) = v(1).*v(2);
A(2) = v(1) .* v(2).^2 + v(4);
A(3) = v(1) .* v(3) .* v(2);
gamma = 1.4;
R = 287; % J/KgK
cv = R/(gamma-1);

e = cv*v(4)./(v(1)*R) + .5*(v(2).^2 + v(3).^2);

A(4) = (e + v(4) ./ v(1)) .* v(1) .* v(2);


