function c = multiply_quatrot(a,b)
a0 = a(:,1); b0 = b(:,1);
a1 = a(:,2); b1 = b(:,2);
a2 = a(:,3); b2 = b(:,3);
a3 = a(:,4); b3 = b(:,4);

c0 = a0.*b0 - a1.*b1 - a2.*b2 - a3.*b3;
c1 = a0.*b1 + a1.*b0 + a2.*b3 - a3.*b2;
c2 = a0.*b2 + a2.*b0 - a1.*b3 + a3.*b1;
c3 = a0.*b3 + a3.*b0 + a1.*b2 - a2.*b1;

c  = [c0,c1,c2,c3];
end
