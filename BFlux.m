function B = BFlux(v)
B = zeros(size(v,1),size(v,2),size(v,3));
B(1) = v(1).*v(3);
B(3) = v(1) .* v(3).^2 + v(4);
B(2) = v(1) .* v(3) .* v(2);
gamma = 1.4;
R = 287; % J/KgK
cv = R/(gamma-1);

e = cv*v(4)./(v(1)*R) + .5*(v(2).^2 + v(3).^2);

B(:,:,4) = (e + v(4) ./ v(1)).*v(1).*v(3);