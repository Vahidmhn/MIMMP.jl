function FindIntersectingPoint(X1,X2,Plane)
    D = length(X1);
    dX = X2-X1;
    Ind = find(dX.!=0);
    Ind = Ind[1];
    t = (-Plane[end]-sum(Plane[1:end-1].*X1))/sum(Plane[1:end-1].*dX);
    X = t.*dX+X1;
end
