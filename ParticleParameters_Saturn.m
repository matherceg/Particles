% Calculate Particle gyro radius, bounce and drift period
% Matija Herceg, December 2021
clear
clc

% ================================
% ================================
% Pitch angle
pitchAngle=90;
% Particle is injected in the magnetic equator at the distance Inj_R
% Testing results from Thomsen and Van Allen, 1980
Inj_R=3.092;
% Injected Energy
nrg=10;   % MeV
% ================================
% ================================



cd('E:\Missions\JUNO\DriftShell')

% Load Mag model
modelName = 'E:\Missions\JUNO\MagneticModel\JupiterGlobalField\Saturn_Z3.mat' ;
load(modelName);
% Max degree
n_max=3;
% Extract model up to the highest degree
model=Saturn_Z3(1:n_max*(n_max+2));

particle=1;

if particle==1
    q=-1.6021766208*10^-19; % Electric charge of electron [C]
    m0=9.10938356*10^-31; % Mass of electron [kg]
elseif particle==2
    q=1.6021766208e-19;    % Electric charge of proton [C]
    m0=1.67262192369e-27;   % Mass of proton [kg]
end

c=299792458;        % Speed of Light [m/s]


% Saturn radius
R=58232e3; %m

% Permeability of vacuum
mu0=4*pi*1e-7; %

% Magnetic Moment
mm0=[model(2)*1e-9,model(3)*1e-9,model(1)*1e-9]*R^3*(4*pi/mu0);
% Total magnetic moment
% sqrt(sum(mm0.^2))


%          E      v/c             gamma
% E_table=[ 10,  0.998693546973503,  20;
%           20,  0.999673546809972,  39;
%           30,  0.999854922852228,  59;
%           40,  0.999918396694404,  78;
%           50,  0.999947774651587,  98;
%           60,  0.999963732686320, 117;
%           70,  0.999973354754878, 137;
%           80,  0.999979599797904, 157;
%           90,  0.999983881356308, 176;
%          100,  0.999986943918602, 196];
%            E              v/c             gamma
% E_table=[   10,	0.998817560667070,	20.5695119784913;
%             20,	0.999689612806676,	40.1390239569826;
%             30,	0.999859742007891,	59.7085359354739;
%             40,	0.999920442452884,	79.2780479139652;
%             50,	0.999948826018144,	98.8475598924565;
%             60,	0.999964342646106,	118.417071870948;
%             70,	0.999973739561062,	137.986583849439;
%             80,	0.999979857936682,	157.556095827930;
%             90,	0.999984062846413,	177.125607806422;
%            100,	0.999987076336511,	196.695119784913];
% Results are checked here: http://hyperphysics.phy-astr.gsu.edu/hbase/Relativ/releng.html

% Energy of the particle in MeV
% E=((E_table(nrg,1))*1e6)'*q;
% v=sqrt(1-(1./(E/(m0*c^2))).^2)*c; % in units of c
% v/c
% Recalculate
% gamma=1/sqrt((1-(v/c)^2))
% E= (gamma * m0 * c^2)/q*1e-6 % MeV

E=(nrg*1e6)'*q;
gamma=abs(E/(m0*c^2))+1;
v=c*sqrt(1-(1/gamma^2));

% Relativistic mass
gamma=1/sqrt((1-(norm(v)/c)^2));
mRes=gamma*m0;

% Rest Energy
E0=m0*c^2;
E0_mev=E0/(1e6*q);
% Kinetic Energy
Ek=m0*c^2*(gamma-1);
Ek_mev=Ek/(1e6*q);
% Relativistic Energy
Erel=gamma*m0*c^2;
Erel_mev=Erel/(1e6*q);

velScale=v/c;

% Integration time
dt=5*1e-7;

% Matrix of rotation from JMR03 to S3
% offset from the rotation axis by 10.31?
% toward system 3 longitude of 196.61?.
mat_Mag_2_S = ea2matr(2,3,2,0, 163.39, 10.31,'deg')';



% Saturn
% Particle bouncing period (Guio et al., 2020)
tempE=Ek_mev*(Ek_mev+2*E0_mev);
Tb = 0.804*(Erel_mev/sqrt(tempE))*Inj_R*(1.25-0.49*sin(pitchAngle*pi/180)-0.04*Inj_R*sin(pitchAngle*pi/180));

% Particle bouncing period (Thomsen and Van Allen, 1980) 
Tb = 0.8006*(Erel_mev/sqrt(tempE))*Inj_R*(1.38-0.32*sin(pitchAngle*pi/180)-0.32*sqrt(sin(pitchAngle*pi/180)));

% Particle averaged drift period
Td = 44.71*(Erel_mev/tempE)*1/Inj_R*(0.44-0.25*sin(pitchAngle*pi/180)+0.12*Inj_R*sin(pitchAngle*pi/180)-7.18e-3*Inj_R^2*sin(pitchAngle*pi/180))^-1;

% Gyroperiod (Thomsen and Van Allen, 1980) 
Tg = 3.495e-6*Inj_R^3*(Erel_mev);

% Gyroradius (Thomsen and Van Allen, 1980)
rG = 1.667e4 *(Inj_R)^3 *sin(pitchAngle*pi/180) * sqrt(Ek_mev*(Ek_mev+2*E0_mev));

% Drift per bounce in degrees
DpB=360/(Td/(Tb/3600));
DpB_km=((2*R*Inj_R*pi/1000)/360)*DpB;

disp(['Injected Energy = ',num2str(nrg),' MeV']);
disp(['Injected Radius = ',num2str(Inj_R),' Rj']);
disp(['Rest Energy = ',num2str(E0)]);
disp(['Kinetic Energy = ',num2str(Ek)]);
disp(['Pitch angle = ',num2str(pitchAngle),'deg']);
disp(['Drift period = ',num2str(abs(Td)),' h']);
disp(['Bouncing period = ',num2str(abs(Tb)),' sec']);
disp(['Drift / Bounce = ',num2str(DpB),' deg (',num2str(DpB_km),'km)']);
disp(['Gyro period = ',num2str(abs(Tg)),' sec']);
disp(['Gyro radius = ',num2str(rG*10^-5),' km']);


% Perpendicular velocity
vPerp = velScale*c*sin(pitchAngle*pi/180);

% Surface mag field [T] of Saturn (Thomsen and Van Allen, 1980) divided by the distance^3
B=0.2e-4/Inj_R^3;

% Radius of gyration
rg=(gamma*m0*vPerp)/(abs(q)*B)/1000

%----------

% Radius of gyration using the Saturns model
cd E:\Missions\JUNO\MagneticModel\JunoMagCalc

% Saturns magnetic field is aligned with Saturn's rotation axis 
% (Connerney et al., 1982)
mat_Mag_2_S = ea2matr(2,3,2, 0,0,0,'deg')';

% Injected particle is in the Magnetic Equator
r_S3  = (mat_Mag_2_S*[Inj_R 0 0]')';



% Convert to spherical Saturn coordinates
[lon,lat,rad] = cart2sph(r_S3(1),r_S3(2),r_S3(3));
Inj_S_sph=[lon*180/pi,lat*180/pi,rad];

% Injected particle in Saturn
r0= Inj_S_sph(3);
lon0=Inj_S_sph(1);
lat0=Inj_S_sph(2);

% Calculate Mag vector from the model (in nT)
[b_r, b_theta, b_phi]= synth_grid_planetary(model,r0,90-lat0,lon0);

% Saturns Z3 (in Gauss)
Br=b_r*10^-5;
Btheta=b_theta*10^-5;
Bphi=b_phi*10^-5;

% Saturns Z3 (in Tesla)
Br=b_r*10^-9;
Btheta=b_theta*10^-9;
Bphi=b_phi*10^-9;

B=sqrt(Br^2+Btheta^2+Bphi^2);

% Rotation of B vector from Local to Z3
matR2=ea2matr_v(3,2,1,lon0,90-lat0,0,'deg');
vB=matR2'*([Btheta, Bphi, Br])';

% Add pitch angle relative to vector B
vPitch=ea2matr_v(1,2,3,pitchAngle,pitchAngle,0,'deg')*vB;

% Perpendicular velocity
vPerp = velScale*c*sin(pitchAngle*pi/180);

% Radius of gyration
rg=(gamma*m0*vPerp)/(abs(q)*B)/1000