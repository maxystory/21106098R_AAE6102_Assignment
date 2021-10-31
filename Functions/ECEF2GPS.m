function GPS = ECEF2GPS(Pos)
    x=Pos(1);
    y=Pos(2);
    z=Pos(3);
    % WGS84 ellipsoid constants:
    a = 6378137;
    e = 8.1819190842622e-2;
    % calculations:
    b   = sqrt(a^2*(1-e^2));
    ep  = sqrt((a^2-b^2)/b^2);
    p   = sqrt(x^2+y^2);
    th  = atan2(a*z,b*p);
    lon = atan2(y,x);
    lat = atan((z+ep^2*b*(sin(th))^3)/(p-e^2*a*(cos(th))^3));
    N   = a/sqrt(1-e^2*(sin(lat))^2);
    alt = p/cos(lat)-N;
    % correct for numerical instability in altitude near exact poles:
    % (after this correction, error is about 2 millimeters, which is about
    % the same as the numerical precision of the overall function)
    k=abs(x)<1 & abs(y)<1;
    alt(k) = abs(z(k))-b;
    GPS=[lat,lon,alt];
return