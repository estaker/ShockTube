function u = CFP(v)
%calculate conservative state vector from primitive state vector
u = zeros(size(v,1),size(v,2),size(v,3));
u(:,:,1) = v(:,:,1);
u(:,:,2) = v(:,:,1) .* v(:,:,2);
u(:,:,3) = v(:,:,1) .* v(:,:,3);
gamma = 1.4;
R = 287; % J/KgK
cv = R/(gamma-1);

e = cv*v(:,:,4)./(v(:,:,1)*R) + .5*(v(:,:,2).^2 + v(:,:,3).^2);

u(:,:,4) = v(:,:,1).*e;