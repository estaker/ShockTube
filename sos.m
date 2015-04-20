function s = sos( v )
gamma = 1.4;
R = 287;
T = v(4)/(v(1)*R);
s = sqrt(gamma*R*T);

end

