function [p1 p2] = FindPointLS0test(Xe,LSe)

np = find (LSe >= 0);
nn = find (LSe < 0);

if length(np) == 1
    Xe  = Xe([np ; nn],:);
    LSe = LSe([np ; nn]);
else
    Xe  = Xe([nn ; np],:);
    LSe = LSe([nn ; np]);
end

A = [Xe LSe];
a = A\[1; 1; 1];

p1 = - [a(1), a(2); Xe(2,2) - Xe(1,2), -Xe(2,1) + Xe(1,1)]\[-1; -(Xe(2,2) - Xe(1,2))*Xe(1,1) + (Xe(2,1) - Xe(1,1))*Xe(1,2)];
p2 = - [a(1), a(2); Xe(3,2) - Xe(1,2), -Xe(3,1) + Xe(1,1)]\[-1; -(Xe(3,2) - Xe(1,2))*Xe(1,1)+ (Xe(3,1) - Xe(1,1))*Xe(1,2)];

end