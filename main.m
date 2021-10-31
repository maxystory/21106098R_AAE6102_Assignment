% AAE6102 Assignment 1
%Lee Max Jwo Lem - 21106098R

clc;
clear;
close all;


%Load Data
load('rcvr.dat');
load('eph.dat');

%Matrix storing satellite ECEF position
Pos_xyz_Mat = zeros(32,13);

%Set constants
Mu = 3.986005*10^14;
Omega_dote = 7.2921151467*10^-5;
C=299792458;
F=-4.442807633*10^(-10);
T_amb=20;
P_amb=101;
P_vap=.86;

%A. Calculate XYZ positions for all valid satellite
for ST_No=1:8
    svid = eph(ST_No,2);
    toc = eph(ST_No,3);
    toe = eph(ST_No,4);
    af0 = eph(ST_No,5); 
    af1 = eph(ST_No,6);
    af2 = eph(ST_No,7);
    ura = eph(ST_No,8);
    e = eph(ST_No,9);
    A = eph(ST_No,10)*eph(ST_No,10);
    dn = eph(ST_No,11);
    M0 = eph(ST_No,12);
    w = eph(ST_No,13);
    omg0 = eph(ST_No,14);
    i0 = eph(ST_No,15);
    odot = eph(ST_No,16);
    idot = eph(ST_No,17);
    cus = eph(ST_No,18);
    cuc = eph(ST_No,19);
    cis = eph(ST_No,20);
    cic = eph(ST_No,21);
    crs = eph(ST_No,22);
    crc = eph(ST_No,23);
    iod = eph(ST_No,24);
    pr = rcvr(find(rcvr(:,2)==svid),3);
    cycles = rcvr(find(rcvr(:,2)==svid),4);
    phase = rcvr(find(rcvr(:,2)==svid),5);
    slp_dtct = rcvr(find(rcvr(:,2)==svid),6);
    snr_dtct = rcvr(find(rcvr(:,2)==svid),7);

    t = eph(ST_No,1) - pr/C; %Transmission Time


    tk = t-toe; %Time from ephemeris epoch
    if tk > 302400  % half week compensation
        tk = tk - 604800;
    end
    if tk < -302400  % half week compensation
        tk = tk + 604800;
    end
    
    %Compute corrected mean motion
    %Input:
        %A(meter): Semi-major axis
        %Mu(meter^3/sec^2): value of Earth's universal gravitational parameters
        %tk(sec):Time from Ephemeris Reference Epoch
        %dn(semi-circle/sec):Mean Motion Diference From Computed Value
        %M0(semi-circles)
    %Output:
        %n0(Rad/sec):Computed mean motion 
        %n(semi-circle/sec):Corrected mean motion
        %Mk(semi-circle):Mean anomaly
    n0=sqrt(Mu/(A^3));
    n=n0+dn;
    Mk=M0+n*tk;


    %Solve Kepler Eq. with Initial value Ek0=Mk
    %Use Fzero Function to solve Kepler Eq.
    %Input:
        %Mk(semicircle): Mean Anomaly
        %e:  Eccentricity
    %Output:
        %Ek(Rad):  Eccentric Anomaly
    Ek=fzero(@(x) Kepler_Eq(x,Mk,e),Mk);

    
    %Compute True Anomaly as a Function of Eccentric Anomaly
    %Input:
        %Ek(Rad):  Eccentric Anomaly
        %e      :  Eccentricity
    %Output:
        %nuk(Rad): True Anomaly
    Sin_nuk=sqrt(1-e^2)*sin(Ek)/(1-e*cos(Ek));
    Cos_nuk=(cos(Ek)-e)/(1-e*cos(Ek));
    nuk=atan2(Sin_nuk,Cos_nuk);


    %Compute Argumet of Lattitude and Correction
    %Input:
        %nuk(Rad): True Anomaly
        %w(semicircle):  Argument of perigee
    %Output:
        %Phik(Rad):  Argument of lattitude
    Phik=nuk+w;




    %Compute Argumet of Lattitude,Radious and Inclination Correction
    %SEcond Harmonic Perturbation
    %Input:
        %Phik(Rad):  Argument of lattitude
        %cus:Amplitude of the cosine harmonic correction term to the argument
        %     of lattitude(Rad)
        %cuc:Amplitude of the sine harmonic correction term to the argument
        %     of lattitude(Rad)
        %crs:Amplitude of the cosine harmonic correction term to orbit
        %radius(meter)
        %crc:Amplitude of the sine harmonic correction term to orbit
        %radius(meter)
        %cis:Amplitude of the cosine harmonic correction term to the angle of
        %     inclination(Rad)
        %cic:Amplitude of the sine harmonic correction term to the angle of
        %     inclination(Rad)
    %Output:
        % Delta_uk(Rad):Argument of lattitude correction
        % Delta_rk(meter):Argument of radius correction
        % Delta_ik(Rad):Argument of inclination correction
    Delta_uk=cus*sin(2*Phik)+cuc*cos(2*Phik);
    Delta_rk=crs*sin(2*Phik)+crc*cos(2*Phik);
    Delta_ik=cis*sin(2*Phik)+cic*cos(2*Phik);



    %Compute Corrected Value of Lattitude,Radious,Inclination and Ascendong
    %node
    %Input:
        % Phik(Rad):  Argument of lattitude
        % A(meter): Semi-major axis
        % Ek(Rad):  Eccentric Anomaly
        % e:  Eccentricity    
        % i0(semi-circle): inlination angle at reference time
        % omg0(semi-circle): Refernce Longitude of Ascending Node
        % odot(semi-circle/sec): rate of right ascension
        % Omega_dote(rad/sec):WGS 84 value of the earth's rotation rate
        % idot(semi-circle/sec): rate of inlination angle
        % Delta_uk(Rad):Argument of lattitude correction
        % Delta_rk(meter):Argument of radius correction
        % Delta_ik(Rad):Argument of inclination correction
    %Output:
        % uk(Rad):Corrcted argument of lattitude
        % rk(meter): corrected radius
        % ik(Rad):corrected inclination
        % Omegak(Rad):Corrected longitude of ascending node.
    uk=Phik+Delta_uk;              %Latitude
    rk=A*(1-e*cos(Ek))+Delta_rk;   %Radious
    ik=i0+Delta_ik+idot*tk;  %Inclination
    Omegak=omg0+(odot-Omega_dote)*tk-Omega_dote*toe;

    


    %Compute satellite vehicle position
    %Satellite position in orbital plane
    x_Perk=rk*cos(uk);
    y_Perk=rk*sin(uk);
    %Satellite Position in ECEF
    xk=x_Perk*cos(Omegak)-y_Perk*cos(ik)*sin(Omegak);
    yk=x_Perk*sin(Omegak)+y_Perk*cos(ik)*cos(Omegak);
    zk=y_Perk*sin(ik);
    %Output:
        % xk,yk,zk(meter): Satellite Position in ECEF
        % Mk(semicircle): Mean Anomaly
        % Ek(Rad):  Eccentric Anomaly
        % nuk(Rad): True Anomaly
        % Phik(Rad):  Argument of lattitude
        % uk(Rad):Corrcted argument of lattitude
        % rk(meter): corrected radius
        % ik(Rad):corrected inclination
        % Omegak(Rad):Corrected longitude of ascending node.
        % A(meter):semi-major orbit axis
        % e: orbit eccentricity
    Pos_xyz=[xk yk zk Mk Ek nuk Phik uk rk ik Omegak A e];

    % earth rotation correction
    % rotation angle
    omegatau = Omega_dote*pr/C;
    % rotation matrix
    R = [ cos(omegatau)    sin(omegatau)   0;
          -sin(omegatau)    cos(omegatau)   0;
           0                0               1];
	xyz_rot = R * [xk, yk, zk]';

    Pos_xyz_Mat(svid,:)=[xyz_rot' Mk Ek nuk Phik uk rk ik Omegak A e]; 
end

%B.Compute Satellite Offset Clock Error
for ST_No = 1:8
    svid = eph(ST_No,2);
    toc = eph(ST_No,3);
    af0 = eph(ST_No,5); 
    af1 = eph(ST_No,6);
    af2 = eph(ST_No,7);
    e = eph(ST_No,9);
    A = eph(ST_No,10)*eph(ST_No,10);
    pr = rcvr(find(rcvr(:,2)==svid),3);

    t = eph(ST_No,1) - pr/C; %Transmission Time

    dtr = F*e*sqrt(A)*sin(Pos_xyz_Mat(svid,5));         %relativistic correction
    dts(ST_No,1) = af0+af1*(t-toc)+af2*(t-toc)^2 + dtr;     %code phase offset by polynomial coefficients



end


%Initial Receiver Position
inipos = [-2694685.473, -4293642.366, 3857878.924];
Xu= -2694685.473;
Yu= -4293642.366;
Zu= 3857878.924;
%Initial Clock offset
dtR= 0;
%Matrix to save estimated receiver position
Pos_Rcv=[Xu Yu Zu dtR];
%Initial delta X
GT = [-2700400, -4292560, 3855270];  % target receiver position


%Start Iteration
Iter = 1;
while 1

    %C. Compute Tropospheric Error Correction (H Model)
    S = size(Pos_xyz_Mat);
    m = S(1);
    n = S(2);
    for i=1:m
      [E,A0]=Calc_Azimuth_Elevation(Pos_Rcv(1:3),Pos_xyz_Mat(i,1:3));
      El(i)=E;                                                   %Elevation Rad
      A(i)=A0;                                                    %Azimoth Rad
    end
    %Zenith Hydrostatic Delay
    Kd=1.55208*10^(-4)*P_amb*(40136+148.72*T_amb)/(T_amb+273.16);

    %Zenith Wet Delay
    Kw=-.282*P_vap/(T_amb+273.16)+8307.2*P_vap/(T_amb+273.16)^2;

    for i=1:m
      Denom1(i)=sin(sqrt(El(i)^2+1.904*10^-3));
      Denom2(i)=sin(sqrt(El(i)^2+.6854*10^-3));
      %Troposhpheric Delay Correctoion
      Delta_R_Trop(i)=Kd/Denom1(i)+Kw/Denom2(i);                        % Meter
    end

    

    %D. Compute H Matrix
    for ST_No = 1:8
        svid = eph(ST_No,2);
        r = norm([Pos_xyz_Mat(svid,1) Pos_xyz_Mat(svid,2) Pos_xyz_Mat(svid,3)]-[Xu Yu Zu]);
        ax = (Xu-Pos_xyz_Mat(svid,1))/r;
        ay = (Yu-Pos_xyz_Mat(svid,2))/r;
        az = (Zu-Pos_xyz_Mat(svid,3))/r;
        H(ST_No,:) = [ax, ay, az, C];
    end

    %E. Compute delta pr
    for ST_No = 1:8
        svid = eph(ST_No,2);
        pr = rcvr(find(rcvr(:,2)==svid),3);  %Measure Pseudorange
        r = norm([Pos_xyz_Mat(svid,1) Pos_xyz_Mat(svid,2) Pos_xyz_Mat(svid,3)]-[Xu Yu Zu]); %Geometric Distance
        b(ST_No) = r + C*(dtR-dts(ST_No,1))+Delta_R_Trop(svid);     %Approx Pseudorange
        deltapr(ST_No) = pr-b(ST_No);
    end


    %F. Computer delta X
    deltax(:,Iter) = (H.'*H)^-1 * H.' * deltapr';


    %G. Residual
    for ST_No = 1:8
        svid = eph(ST_No,2);
        pr = rcvr(find(rcvr(:,2)==svid),3);  %Measure Pseudorange
        prmatrix(ST_No,1) = pr;
    end
    residual = prmatrix-(b'+H*deltax(:,Iter)); % calculate residual
    residualSE = residual' * residual; % residual, squared error

    %Update solution
    Xu= Xu + deltax(1,Iter);
    Yu= Yu + deltax(2,Iter);
    Zu= Zu + deltax(3,Iter);
    dtR= dtR + deltax(4,Iter);
    Pos_Rcv(Iter,:)=[Xu Yu Zu dtR];

    posErr = norm([Xu Yu Zu]-GT); % positioning error to given target position

    if norm(deltax(:,Iter)) < 1
        break
    end
    Iter = Iter+1;
end


figure(1)
[ilat, ilon, ialt] = Wgsxyz2lla(inipos); % initial position to wgs84 LLA
[wlat, wlon, walt] = Wgsxyz2lla(Pos_Rcv(Iter,1:3)); % convert updated solution to wgs84 LLA
[glat, glon, galt] = Wgsxyz2lla(GT); % convert GT to wgs84 LLA
 geoscatter([ilat],[ilat],'filled','MarkerFaceColor','black')
 hold on
 geoscatter([glat],[glon],'filled','MarkerFaceColor','red')
 geoscatter([wlat],[wlon],'filled','MarkerFaceColor','yellow')
 legend('Initial Position','Ground Truth','Estimated Position after 5 Iterations')

figure(2)
scatter([ilat],[ilat],500,'filled','MarkerFaceColor','black')
hold on
scatter([glat],[glon],500,'filled','MarkerFaceColor','red')
scatter([wlat],[wlon],100,'filled','MarkerFaceColor','yellow')
legend('Initial Position','Ground Truth','Estimated Position after 5 Iterations')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)') 