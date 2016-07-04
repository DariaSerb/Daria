function Epsilon_matter = AreaRatio(Elements,Number_of_CutElt,Nodes,LS)

Te = Elements(Number_of_CutElt,2:end-2);
LSe = LS(Te);
coord = Nodes(Te,2:end-1);
FullArea = TriArea(coord);

np = find(LSe > 0);
%  [p1 p2] = FindPointLS0test(coord,LSe);
if length(np) == 1
    nn = setdiff(1:3, np);
    % C = setdiff(A,B) returns the data in A that is not in B.
    p1 = FindPointLS0(coord(nn(1),:),LSe(nn(1)),coord(np,:),LSe(np));
    p2 = FindPointLS0(coord(nn(2),:),LSe(nn(2)),coord(np,:),LSe(np));
    
    ElTriArea = TriArea([p1;p2;coord(1,:)]);
    
else
    nn = setdiff(1:3, np);
    p1 = FindPointLS0(coord(np(1),:),LSe(np(1)),coord(nn,:),LSe(nn));
    p2 = FindPointLS0(coord(np(2),:),LSe(np(2)),coord(nn,:),LSe(nn));
    
    ElTriArea = TriArea([p1;p2;coord(2,:)]);
    ElTriArea = FullArea - ElTriArea;
end
Epsilon_matter = ElTriArea/FullArea;
end


function p0 = FindPointLS0(p1,d1,p2,d2)

p0 = p1 + (p2 - p1)*(abs(d1)/(abs(d1) + abs(d2)));

end

function Area = TriArea(X)

Area = 0.5*abs((X(1,1) - X(3,1))*(X(2,2) - X(3,2)) - (X(2,1) - X(3,1))*(X(1,2) - X(3,2)));

end

