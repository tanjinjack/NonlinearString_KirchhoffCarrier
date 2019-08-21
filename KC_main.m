%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the "main" routine for computing displacement
% via Kirchhoff-Carrier equations, expressed in modal form,
% i.e. U = \sum_n^N u_n(t) \phi_n(x). 
% The temporal component, u_n(t) is solved via an energy-conserved
% finite difference sheme and the mode shape, \phi_n(x) is
% simply sin (n \pi x /L)
% refer to J.J. Tan's thesis (2017) Section 2.4.2 (pg 33) for more info
% open access link: https://pastel.archives-ouvertes.fr/tel-01755039
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

CH = load('stringcharacterA22-16.mat');

%%physical parameter
phyParam.rho       = 7850; %density
phyParam.diameter  = 1.3e-3; %string diameter
phyParam.young     = 1.897e11; %Young's modulus
phyParam.T0        = 895.2932; %tension
phyParam.L         = 0.668; %length
phyParam.x0        = 0.334; %absolute plucking location in meter
phyParam.amplitude = 7e-3; %plucking amplitude in meter
phyParam.angle     = 2; %plucking angle (deviation from vertical plane)
phyParam.delta_f   = 0.15; %detuning between two polarisations
phyParam.B         = CH.B; % inharmonicity

phyParam.Qte       = CH.Qte10; %Valette & Cuesta damping parameter
phyParam.deltave   = CH.deltave10; %Valette & Cuesta damping parameter
phyParam.RM        = CH.RM10; % Bensa's damping parameter
phyParam.etaM      = CH.etaM10; % Bensa's damping parameter

%%simulation parameters
simParam.HALFSINE     = 0; %halfsine IC if =1
simParam.Nmode        = 40; %number of modes to compute
simParam.alpha        = 0.5; %newmark scheme parameter
simParam.Tend         = 10; %end time for simulation
simParam.timestep     = 1e-5; % time step for simulation
simParam.NInitCond    = simParam.Nmode; %number of modes for init. cond.
simParam.ObservePoint = 0.638; %absolute location in meter

%%optional parameters
optParam.linearOnly = 0; %0 for nonlinear, 1 for linear

%useDampingModel: 0 for all V&C, 1 for experimental data + V&C, 
%2 for Bensa, 3 for experimental data + Bensa, -1 for NO damping at all
optParam.useDampingModel =  1; 
optParam.dampExpData = -CH.mean_amo; %only be used if useDampingModel == 1 or 3


OUT = KirchhoffCarrier(phyParam,simParam,optParam);
save('KC_comparewithmultamp_VCdamp.mat')

return
