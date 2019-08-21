function OUT = KirchhoffCarrier(phyParam,simParam,optParam)

if(nargin < 3)
    linearOnly = 0;
    useDampingModel = -1;
else
    linearOnly = optParam.linearOnly;
    useDampingModel = optParam.useDampingModel;
end

if useDampingModel == 1 || useDampingModel == 1
    

    if isfield(optParam,'dampExpData') == 0
        disp('You have selected to use experimental damping data.')
        disp('However you did not provide any experimental data.')
        disp('Please make sure the data are in optParam.dampExpData.')
        disp('Exiting code...')
        
        OUT = 0;
        return
    end
end

%preset all damping parameters to 0
%if there is no input to these parameters, they will be 0 by default
Qte      = 0;
deltave  = 0;
RM       = 0;
etaM     = 0;


%%assigning physical parameters
rho      = phyParam.rho;
diameter = phyParam.diameter;
area     = diameter^2*pi/4;
young    = phyParam.young;
I        = pi*diameter^4/64;
T0       = phyParam.T0;
L        = phyParam.L;
x0       = phyParam.x0; %absolute location in meter
amplitude= phyParam.amplitude;
angle    = phyParam.angle;
delta_f  = phyParam.delta_f; 
B        = phyParam.B;
c        = sqrt(T0/rho/area);

Qte      = phyParam.Qte;
deltave  = phyParam.deltave;
RM       = phyParam.RM;
etaM     = phyParam.etaM;

%%assigning simulation parameters
Nmode = simParam.Nmode;
alpha = simParam.alpha;
Tend     = simParam.Tend;
timestep = simParam.timestep;
NInitCond = simParam.NInitCond;
ObservePoint = simParam.ObservePoint;
HALFSINE = simParam.HALFSINE;

%%logic check before proceeding
if ObservePoint > L 
    disp('ObservePoint should be smaller than L, and is expressed in metre.');
    disp('Please check your input. Exiting code now...')
    OUT = 0;
    return
elseif ObservePoint < 0
    disp('ObservePoint should be larger than 0.')
    disp('Please check your input. Exiting code now...')
    OUT = 0;
    return    
end
if x0 > L 
    disp('x0 should be smaller than L, and is expressed in metre.');
    disp('Please check your input. Exiting code now...')
    OUT = 0;
    return
elseif x0 < 0
    disp('x0 should be larger than 0.')
    disp('Please check your input. Exiting code now...')
    OUT = 0;
    return      
end
if Nmode < NInitCond
    disp('There are more modes for NInitCond than for Nmode.')
    disp('NInitCond should always be smaller or equal to Nmode')
    disp('Please check your input. Exiting code now...')
    OUT = 0;
    return
end


Ntimestep = round(Tend/timestep);

spacesteptemp = c*timestep;
Nelement = floor(L/spacesteptemp);

spacestep = L/Nelement;
lambda = c*timestep/spacestep;



%%precalculated variable

n=1:Nmode;

T0v = rho * area * (2*L*delta_f + sqrt(T0/rho/area))^2;
f0 = 1/2/L * sqrt(T0/rho/area);
f0v = 1/2/L * sqrt(T0v/rho/area);
freq = (n*f0+n.^3*B);
freqv = (n*f0v+n.^3*B);
% OHM2 = (freq*2*pi).^2;
% OHM2v = (freqv*2*pi).^2;

OHM2 = young*I*n.^4*pi^4/(rho*area*L^4) + T0*n.^2*pi^2/(rho*area*L^2);
OHM2v = young*I*n.^4*pi^4/(rho*area*L^4) + T0v*n.^2*pi^2/(rho*area*L^2);
if(freq(end) > (1/timestep)/2)
    disp('WARNING:')
    disp(['The ',num2str(Nmode),'-th frequency is ',num2str(freq(end)),'Hz '])
    disp(['but sampling frequency is only up until ',num2str((1/timestep)),'Hz '])
    disp('Sampling frequency needs to be 2 times more than the highest frequency')
    disp('Potential issue in computation')
end
%%damping%%
etaair = 1.8e-5;
rhoair = 12;
Qair_m1 = (2*pi*etaair+2*pi*diameter.*sqrt(pi*etaair*rhoair*freq)) ./ (2*pi*rho*area * freq);
Qve_m1  = 4*pi^2*rho*area*young*I*deltave/T0^2 * freq.^2;
Qte_m1  = 1/Qte;
sigma = pi*freq.*(Qair_m1+Qve_m1+Qte_m1);

if useDampingModel == 1
    whicheversmaller = min([length(optParam.dampExpData) Nmode]);
    sigma(1:whicheversmaller) = optParam.dampExpData(1:whicheversmaller);
elseif useDampingModel == 2
    sigma = RM + etaM*(2*pi*freq).^2;
elseif useDampingModel == 3
    sigma = RM + etaM*(2*pi*freq).^2;
    whicheversmaller = min([length(optParam.dampExpData) Nmode]);
    sigma(1:whicheversmaller) = optParam.dampExpData(1:whicheversmaller);
elseif useDampingModel == -1
    sigma = zeros(size(sigma));
end

Q = 2*sigma;

%nonlinear parameters
if linearOnly == 1
    N = 0;
else
    N = young * n.^2 * pi^4 / 2 / rho / L^5 *L/2 /2;
end


A = 1/timestep^2 + OHM2*(1-alpha)/2 + Q/2/timestep;
D = 1/timestep^2 + OHM2v*(1-alpha)/2 + Q/2/timestep;


%%initial condition
indexvec = 1:NInitCond;
IC = 2/L./(indexvec*pi/L).^2 .* ...
    (1/x0*(sin(indexvec*pi*x0/L)-indexvec*pi*x0/L.*cos(indexvec*pi*x0/L))+ ...
    1/(L-x0)*(indexvec*pi/L*(L-x0).*cos(indexvec*pi*x0/L)+sin(indexvec*pi*x0/L)));

ICu = IC;
ICv = IC;
% IC
if HALFSINE == 1
    IC = IC*1e-16;
    IC(1) = 1;
    ICu = IC;
    ICv = 1e-16*IC;
end
    

u0 = zeros([1 Nmode]);
v0 = zeros([1 Nmode]);
u0(1:NInitCond) = amplitude*cos(angle*pi/180)*ICu;
v0(1:NInitCond) = amplitude*sin(angle*pi/180)*ICv;

%%set initial condition that is 0 to 1e-16
u0zero = find(u0==0);
for i = 1:length(u0zero)
    u0(u0zero(i)) = 1e-16;
end
v0zero = find(v0==0);
for i = 1:length(v0zero)
    v0(v0zero(i)) = 1e-16;
end
 
%%initialise u and v
U_np1 = zeros([Ntimestep Nmode]);
U_n   = zeros([Ntimestep Nmode]);
U_nm1 = zeros([Ntimestep Nmode]);

V_np1 = U_np1;
V_n   = U_n;
V_nm1 = U_nm1;

U_np1(1,:) = u0;
U_n(1,:)   = u0;
U_nm1(1,:) = u0;

V_np1(1,:) = v0;
V_n(1,:)   = v0;
V_nm1(1,:) = v0;

%%initialise J and K matrix
J = zeros(Nmode*2);
K = zeros([Nmode 1]);
uv_sum = zeros([1 Ntimestep]);
ENERGY = zeros([1 Ntimestep]);


%%the loop
for i = 2:Ntimestep
    
    U_n(i,:) = U_np1(i-1,:);
    V_n(i,:) = V_np1(i-1,:);
    U_nm1(i,:) = U_n(i-1,:);
    V_nm1(i,:) = V_n(i-1,:);
    

    B = (N.*U_n(i,:))'*((n.^2).*U_n(i,:));
    C = (N.*U_n(i,:))'*((n.^2).*V_n(i,:));
    E = (N.*V_n(i,:))'*((n.^2).*U_n(i,:));
    F = (N.*V_n(i,:))'*((n.^2).*V_n(i,:));

    J = [[B+diag(A) C]; [E F+diag(D)];];

    
    uv_sum(i) = sum(n.^2.*(U_nm1(i,:).*U_n(i,:) + V_nm1(i,:).*V_n(i,:)));
    G = -((-2/timestep^2 + alpha*OHM2+N*uv_sum(i)).*U_n(i,:) ...
        + (1/timestep^2 + OHM2*(1-alpha)/2 - Q/2/timestep).*U_nm1(i,:));
    H = -((-2/timestep^2 + alpha*OHM2v+N*uv_sum(i)).*V_n(i,:) ...
        + (1/timestep^2 + OHM2v*(1-alpha)/2 - Q/2/timestep).*V_nm1(i,:));    
    K = [G H]';
    
    W = J\K;
    
    U_np1(i,:) = W(1:Nmode);
    V_np1(i,:) = W((Nmode+1):end);
    
    NRG1 = ((U_n(i,:) - U_nm1(i,:))/timestep).^2;
    NRG2 = ((V_n(i,:) - V_nm1(i,:))/timestep).^2;
    
    NRG3 = ((U_n(i,:) + U_nm1(i,:))/2).^2;
    NRG4 = ((V_n(i,:) + V_nm1(i,:))/2).^2;
    
    NRGNL = uv_sum(i).^2;
    
    coef1 = (1+(1-2*alpha)/4*timestep^2*OHM2);
    coef2 = (1+(1-2*alpha)/4*timestep^2*OHM2v);
    coef3 = OHM2;
    coef4 = OHM2v;
    
    sumNRG1 = sum(coef1.*NRG1);
    sumNRG2 = sum(coef2.*NRG2);
    sumNRG3 = sum(coef3.*NRG3);
    sumNRG4 = sum(coef4.*NRG4);

        
    
    if linearOnly == 1
        coefNL = 0;
    else
        coefNL = young*pi^4/8/rho/L^4;
    end
    
      
    ENERGY(i) = sumNRG1 + sumNRG2 + sumNRG3 + sumNRG4 + coefNL*NRGNL;
       
end
x = ObservePoint;
phi =  sin(n*pi*x/L);

U = U_n*phi';
V = V_n*phi';
theta = atan2(V,U);
Tvec = 0:timestep:Tend-timestep*1;

dE = ENERGY(4:end)-ENERGY(3:end-1);
% figure;plot(dE);title('dE')
% figure;plot(ENERGY)
OUT.T = Tvec;
OUT.U = U;
OUT.V = V;
OUT.NRG = ENERGY;

end