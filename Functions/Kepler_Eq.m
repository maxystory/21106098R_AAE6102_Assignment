%Kepler Eq => M=E-e*sinE

function f=Kepler_Eq(x,M,e)
f=x-e*sin(x)-M;