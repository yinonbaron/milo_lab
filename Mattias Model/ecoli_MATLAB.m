
function sim = ecoliOK

%% Description of the file
%
% This file contains a differential equation model of E. coli's
% central carbon metabolism and its regulation.
% The model is capable to adapt on the mechanistic scale from
% growth on glucose to growth on acetate and vice versa.
% The model is presented in detail in the paper:
%
% 'Bacterial adaptation through distributed sensing of metabolic fluxes'
% by Oliver Kotte, Judith B. Zaugg & and Matthias Heinemann,
% In press, Molecular Systems Bioloy, 2010
% 
% To run the model, type
% sim = ecoliOK;
% into the MATLAB Command Window. The function ecoliOK then returns
% the simulation output as a structure 'sim', with
% sim.t the time points of the simulation
% sim.x the simulated states
% sim.f the simulated rates
% sim.p the used parameters
% sim.ic the used initial conditions
% sim.S the stoichiometric matrix
% sim.a a structure of aliases to easily access values of interest
% sim.a the chosen simulation scenario (see below)
%
% The returned structure sim can be used to plot and analyze the simulation
% results. For instance, the simulated Cra-FBP levels over time
% are plotted by typing
% plot ( sim.t, sim.x(:,sim.a.x.CraFBP) )
% into the MATLAB Command Window.
% For a definition of all aliases, see below.
%
% To specify the simulation, choose one of the six predefined scenarios
% of interest, or define a custom simulation (see below).
%
% Internally, the function ecoliOK follows this sequential workflow:
% The stoichiometric matrix S is defined.
% The parameter values p are defined.
% The initial conditions and the simulation time are defined.
% The solver ode15s calls the nested function ecoliOK_core, which
%   first computes the rates f as a function of p and x (the states),
%   and then uses these to compute the time derivatives dx of the
%   states, through dx = S * f.
% The simulation output is processed and stored in the structure sim,
%   which forms the output of the main function ecoliOK.


%% Choose the scenario to be simulated
% To run a simulation, choose one of the predefined scenarios 1-6,
% or define a custom simulation. To define a custom simulation,
% set 'scenario=7', search for the string 'if scenario==7,'
% and enter the desired simulation time and initial conditions there.
% To customize the parameter values, search for the string
% 'Definition of the parameter values' and change the values there.
%
% The available scenarios are:
% 1: The population is adapted to glucose, the carbon source is glucose
% 2: The population is adapted to acetate, the carbon source is acetate
% 3: The population is adapted to glucose, the carbon source is acetate
% 4: The population is adapted to acetate, the carbon source is glucose
% 5: Diauxic shift scenario (the first environment of scenario 6)
% 6: Extended diauxic shift scenario (3 consecutive environments),
%    as presented in the paper
% 7: Custom simulation

scenario = 4;


%% Initialize variables

a = [];
S = [];
p = [];
ssACT = [];
ssGLC = [];
ic = [];
sim = [];


%% Definition of Aliases
% The purpose of the aliases is to make the MATLAB code readable,
% and to make the simulation output accessible.
% All aliases are stored in the structure 'a':
% 'a.x.NAME' denote state variables
% 'a.f.(...)' denote rates
%   'a.f.e.NAME'   -> enzymatic reaction rates
%   'a.f.pts.NAME' -> rates in the PTS system
%   'a.f.tf.NAME'  -> metabolite-transcription factor binding rates
%   'a.f.g.NAME'   -> gene expression rates
%   'a.f.d.NAME'   -> dilution and degradation rates
%   'a.f.bm.NAME'  -> biomass-producing reaction rates
%   'a.f.env.NAME' -> interaction rates with the cell's environment
% 'a.p.(...)' denote parameters
%   'a.p.e.NAME    -> enzyme kinetics
%   'a.p.pts.NAME  -> PTS system
%   'a.p.tf.NAME   -> metabolite-transcription factor interactions
%   'a.p.g.NAME    -> gene expression
%   'a.p.d.NAME    -> protein degradation
%   'a.p.bm.NAME   -> growth and biomass production
%   'a.p.env.NAME  -> interaction with the environment

% Definition of state aliases

a.x.OD = 1;
a.x.ACT = 2;
a.x.GLC = 3;
a.x.ACoA = 4;
a.x.AKG = 5;
a.x.cAMP = 6;
a.x.FBP = 7;
a.x.G6P = 8;
a.x.GLX = 9;
a.x.ICT = 10;
a.x.MAL = 11;
a.x.OAA = 12;
a.x.PEP = 13;
a.x.PG3 = 14;
a.x.PYR = 15;
a.x.AceA = 16;
a.x.AceB = 17;
a.x.AceK = 18;
a.x.Acoa2act = 19;
a.x.Acs = 20;
a.x.Akg2mal = 21;
a.x.CAMPdegr = 22;
a.x.Cya = 23;
a.x.Emp = 24;
a.x.Eno = 25;
a.x.Fdp = 26;
a.x.GltA = 27;
a.x.Icd = 28;
a.x.Icd_P = 29;
a.x.Mdh = 30;
a.x.MaeAB = 31;
a.x.PckA = 32;
a.x.Pdh = 33;
a.x.PfkA = 34;
a.x.Ppc = 35;
a.x.PpsA = 36;
a.x.PykF = 37;
a.x.EIIA = 38;
a.x.EIIA_P = 39;
a.x.EIICB = 40;
a.x.Cra = 41;
a.x.CraFBP = 42;
a.x.Crp = 43;
a.x.CrpcAMP = 44;
a.x.IclR = 45;
a.x.PdhR = 46;
a.x.PdhRPYR = 47;

% Definition of rate aliases

a.f.env.growth = 1;
a.f.env.GLCup = 2;
a.f.env.ACTup = 3;
a.f.env.ACTex = 4;
a.f.e.AceA = 5;
a.f.e.AceB = 6;
a.f.e.AceK_Ki = 7;
a.f.e.AceK_Ph = 8;
a.f.e.Acoa2act = 9;
a.f.e.Acs = 10;
a.f.e.Akg2mal = 11;
a.f.e.CAMPdegr = 12;
a.f.e.Cya = 13;
a.f.e.Emp = 14;
a.f.e.Eno = 15;
a.f.e.Fdp = 16;
a.f.e.GltA = 17;
a.f.e.Icd = 18;
a.f.e.MaeAB = 19;
a.f.e.Mdh = 20;
a.f.e.PckA = 21;
a.f.e.Pdh = 22;
a.f.e.PfkA = 23;
a.f.e.Ppc = 24;
a.f.e.PpsA = 25;
a.f.e.PykF = 26;
a.f.pts.r1 = 27;
a.f.pts.r4 = 28;
a.f.tf.Cra = 29;
a.f.tf.Crp = 30;
a.f.tf.PdhR = 31;
a.f.g.aceA = 32;
a.f.g.aceB = 33;
a.f.g.aceK = 34;
a.f.g.acoa2act = 35;
a.f.g.acs = 36;
a.f.g.akg2mal = 37;
a.f.g.campdegr = 38;
a.f.g.cra = 39;
a.f.g.crp = 40;
a.f.g.cya = 41;
a.f.g.EIIA = 42;
a.f.g.EIICB = 43;
a.f.g.emp = 44;
a.f.g.eno = 45;
a.f.g.fdp = 46;
a.f.g.gltA = 47;
a.f.g.icd = 48;
a.f.g.iclr = 49;
a.f.g.mdh = 50;
a.f.g.maeAB = 51;
a.f.g.pckA = 52;
a.f.g.pdh = 53;
a.f.g.pdhr = 54;
a.f.g.pfkA = 55;
a.f.g.ppc = 56;
a.f.g.ppsA = 57;
a.f.g.pykF = 58;
a.f.d.ACoA = 59;
a.f.d.AKG = 60;
a.f.d.cAMP = 61;
a.f.d.FBP = 62;
a.f.d.G6P = 63;
a.f.d.GLX = 64;
a.f.d.ICT = 65;
a.f.d.MAL = 66;
a.f.d.OAA = 67;
a.f.d.PEP = 68;
a.f.d.PG3 = 69;
a.f.d.PYR = 70;
a.f.d.AceA = 71;
a.f.d.AceB = 72;
a.f.d.AceK = 73;
a.f.d.Acoa2act = 74;
a.f.d.Acs = 75;
a.f.d.Akg2mal = 76;
a.f.d.CAMPdegr = 77;
a.f.d.Cra = 78;
a.f.d.CraFBP = 79;
a.f.d.Crp = 80;
a.f.d.CrpcAMP = 81;
a.f.d.Cya = 82;
a.f.d.EIIA = 83;
a.f.d.EIIA_P = 84;
a.f.d.EIICB = 85;
a.f.d.Emp = 86;
a.f.d.Eno = 87;
a.f.d.Fdp = 88;
a.f.d.GltA = 89;
a.f.d.Icd = 90;
a.f.d.Icd_P = 91;
a.f.d.IclR = 92;
a.f.d.Mdh = 93;
a.f.d.MaeAB = 94;
a.f.d.PckA = 95;
a.f.d.Pdh = 96;
a.f.d.PdhR = 97;
a.f.d.PdhRPYR = 98;
a.f.d.PfkA = 99;
a.f.d.Ppc = 100;
a.f.d.PpsA = 101;
a.f.d.PykF = 102;
a.f.bm.ACoA = 103;
a.f.bm.AKG = 104;
a.f.bm.G6P = 105;
a.f.bm.OAA = 106;
a.f.bm.PEP = 107;
a.f.bm.PG3 = 108;
a.f.bm.PYR = 109;

% Definition of parameter aliases

a.p.env.M_ACT = 1;
a.p.env.M_GLC = 2;
a.p.env.uc = 3;
a.p.e.AceA.kcat = 4;
a.p.e.AceA.n = 5;
a.p.e.AceA.L = 6;
a.p.e.AceA.Kict = 7;
a.p.e.AceA.Kpep = 8;
a.p.e.AceA.Kpg3 = 9;
a.p.e.AceA.Kakg = 10;
a.p.e.AceB.kcat = 11;
a.p.e.AceB.Kglx = 12;
a.p.e.AceB.Kacoa = 13;
a.p.e.AceB.Kglxacoa = 14;
a.p.e.AceK.kcat_ki = 15;
a.p.e.AceK.kcat_ph = 16;
a.p.e.AceK.n = 17;
a.p.e.AceK.L = 18;
a.p.e.AceK.Kicd = 19;
a.p.e.AceK.Kicd_P = 20;
a.p.e.AceK.Kpep = 21;
a.p.e.AceK.Kpyr = 22;
a.p.e.AceK.Koaa = 23;
a.p.e.AceK.Kglx = 24;
a.p.e.AceK.Kakg = 25;
a.p.e.AceK.Kpg3 = 26;
a.p.e.AceK.Kict = 27;
a.p.e.Acoa2act.kcat = 28;
a.p.e.Acoa2act.n = 29;
a.p.e.Acoa2act.L = 30;
a.p.e.Acoa2act.Kacoa = 31;
a.p.e.Acoa2act.Kpyr = 32;
a.p.e.Acs.kcat = 33;
a.p.e.Acs.Kact = 34;
a.p.e.Akg2mal.kcat = 35;
a.p.e.Akg2mal.Kakg = 36;
a.p.e.CAMPdegr.kcat = 37;
a.p.e.CAMPdegr.KcAMP = 38;
a.p.e.Cya.kcat = 39;
a.p.e.Cya.KEIIA = 40;
a.p.e.Emp.kcat.f = 41;
a.p.e.Emp.kcat.r = 42;
a.p.e.Emp.Kfbp = 43;
a.p.e.Emp.Kpg3 = 44;
a.p.e.Eno.kcatf = 45;
a.p.e.Eno.kcatr = 46;
a.p.e.Eno.Kpg3 = 47;
a.p.e.Eno.Kpep = 48;
a.p.e.Fdp.kcat = 49;
a.p.e.Fdp.n = 50;
a.p.e.Fdp.L = 51;
a.p.e.Fdp.Kfbp = 52;
a.p.e.Fdp.Kpep = 53;
a.p.e.GltA.kcat = 54;
a.p.e.GltA.Koaa = 55;
a.p.e.GltA.Kacoa = 56;
a.p.e.GltA.Koaaacoa = 57;
a.p.e.GltA.Kakg = 58;
a.p.e.Icd.kcat = 59;
a.p.e.Icd.n = 60;
a.p.e.Icd.L = 61;
a.p.e.Icd.Kict = 62;
a.p.e.Icd.Kpep = 63;
a.p.e.Mdh.kcat = 64;
a.p.e.Mdh.n = 65;
a.p.e.Mdh.Kmal = 66;
a.p.e.MaeAB.kcat = 67;
a.p.e.MaeAB.n = 68;
a.p.e.MaeAB.L = 69;
a.p.e.MaeAB.Kmal = 70;
a.p.e.MaeAB.Kacoa = 71;
a.p.e.MaeAB.Kcamp = 72;
a.p.e.PckA.kcat = 73;
a.p.e.PckA.Koaa = 74;
a.p.e.PckA.Kpep = 75;
a.p.e.Pdh.kcat = 76;
a.p.e.Pdh.n = 77;
a.p.e.Pdh.L = 78;
a.p.e.Pdh.Kpyr = 79;
a.p.e.Pdh.KpyrI = 80;
a.p.e.Pdh.Kglx = 81;
a.p.e.PfkA.kcat = 82;
a.p.e.PfkA.n = 83;
a.p.e.PfkA.L = 84;
a.p.e.PfkA.Kg6p = 85;
a.p.e.PfkA.Kpep = 86;
a.p.e.Ppc.kcat = 87;
a.p.e.Ppc.n = 88;
a.p.e.Ppc.L = 89;
a.p.e.Ppc.Kpep = 90;
a.p.e.Ppc.Kfbp = 91;
a.p.e.PpsA.kcat = 92;
a.p.e.PpsA.n = 93;
a.p.e.PpsA.L = 94;
a.p.e.PpsA.Kpyr = 95;
a.p.e.PpsA.Kpep = 96;
a.p.e.PykF.kcat = 97;
a.p.e.PykF.n = 98;
a.p.e.PykF.L = 99;
a.p.e.PykF.Kpep = 100;
a.p.e.PykF.Kfbp = 101;
a.p.pts.k1 = 102;
a.p.pts.km1 = 103;
a.p.pts.k4 = 104;
a.p.pts.KEIIA = 105;
a.p.pts.Kglc = 106;
a.p.tf.Cra.scale = 107;
a.p.tf.Cra.kfbp = 108;
a.p.tf.Cra.n = 109;
a.p.tf.Crp.scale = 110;
a.p.tf.Crp.kcamp = 111;
a.p.tf.Crp.n = 112;
a.p.tf.PdhR.scale = 113;
a.p.tf.PdhR.kpyr = 114;
a.p.tf.PdhR.n = 115;
a.p.g.aceBAK.vcra_unbound = 116;
a.p.g.aceBAK.vcra_bound = 117;
a.p.g.aceBAK.Kcra = 118;
a.p.g.aceBAK.aceBfactor = 119;
a.p.g.aceBAK.aceKfactor = 120;
a.p.g.aceBAK.KDNA = 121;
a.p.g.aceBAK.KP = 122;
a.p.g.aceBAK.KPprime = 123;
a.p.g.aceBAK.KG = 124;
a.p.g.aceBAK.L = 125;
a.p.g.aceBAK.kcat_iclr = 126;
a.p.g.aceBAK.DNA = 127;
a.p.g.aceBAK.vcrp_bound = 128;
a.p.g.aceBAK.vcrp_unbound = 129;
a.p.g.aceBAK.Kcrp = 130;
a.p.g.acs.vcrp_unbound = 131;
a.p.g.acs.vcrp_bound = 132;
a.p.g.acs.n = 133;
a.p.g.acs.Kcrp = 134;
a.p.g.akg2mal.vcrp_unbound = 135;
a.p.g.akg2mal.vcrp_bound = 136;
a.p.g.akg2mal.Kcrp = 137;
a.p.g.akg2mal.n = 138;
a.p.g.emp.vcra_unbound = 139;
a.p.g.emp.vcra_bound = 140;
a.p.g.emp.Kcra = 141;
a.p.g.emp.vcrp_unbound = 142;
a.p.g.emp.vcrp_bound = 143;
a.p.g.emp.Kcrp = 144;
a.p.g.eno.vcra_unbound = 145;
a.p.g.eno.vcra_bound = 146;
a.p.g.eno.Kcra = 147;
a.p.g.fdp.vcra_unbound = 148;
a.p.g.fdp.vcra_bound = 149;
a.p.g.fdp.Kcra = 150;
a.p.g.gltA.vcrp_unbound = 151;
a.p.g.gltA.vcrp_bound = 152;
a.p.g.gltA.Kcrp = 153;
a.p.g.gltA.n = 154;
a.p.g.icd.vcra_unbound = 155;
a.p.g.icd.vcra_bound = 156;
a.p.g.icd.Kcra = 157;
a.p.g.mdh.vcrp_unbound = 158;
a.p.g.mdh.vcrp_bound = 159;
a.p.g.mdh.Kcrp = 160;
a.p.g.pckA.vcra_unbound = 161;
a.p.g.pckA.vcra_bound = 162;
a.p.g.pckA.Kcra = 163;
a.p.g.pdh.vpdhr_unbound = 164;
a.p.g.pdh.vpdhr_bound = 165;
a.p.g.pdh.Kpdhr = 166;
a.p.g.pfkA.vcra_unbound = 167;
a.p.g.pfkA.vcra_bound = 168;
a.p.g.pfkA.Kcra = 169;
a.p.g.ppsA.vcra_unbound = 170;
a.p.g.ppsA.vcra_bound = 171;
a.p.g.ppsA.Kcra = 172;
a.p.g.pykF.vcra_unbound = 173;
a.p.g.pykF.vcra_bound = 174;
a.p.g.pykF.Kcra = 175;
a.p.d.k_degr = 176;
a.p.bm.k_expr = 177;
a.p.bm.muACT = 178;
a.p.bm.muGLC = 179;
a.p.bm.GLC.ACoA = 180;
a.p.bm.GLC.AKG = 181;
a.p.bm.GLC.G6P = 182;
a.p.bm.GLC.OAA = 183;
a.p.bm.GLC.PEP = 184;
a.p.bm.GLC.PG3 = 185;
a.p.bm.GLC.PYR = 186;
a.p.bm.ACT.ACoA = 187;
a.p.bm.ACT.AKG = 188;
a.p.bm.ACT.G6P = 189;
a.p.bm.ACT.OAA = 190;
a.p.bm.ACT.PEP = 191;
a.p.bm.ACT.PG3 = 192;
a.p.bm.ACT.PYR = 193;

%% Definition of the stoichiometric matrix 'S'

nr_of_states = 47;
nr_of_rates = 109;

S = sparse(nr_of_states,nr_of_rates);

S(a.x.OD,a.f.env.growth) = 1;
S(a.x.GLC,a.f.env.GLCup) = -1;
S(a.x.ACT,a.f.env.ACTup) = -1;
S(a.x.ACT,a.f.env.ACTex) = 1;
S(a.x.AceA,a.f.g.aceA) = 1;
S(a.x.AceA,a.f.d.AceA) = -1;
S(a.x.AceB,a.f.g.aceB) = 1;
S(a.x.AceB,a.f.d.AceB) = -1;
S(a.x.AceK,a.f.g.aceK) = 1;
S(a.x.AceK,a.f.d.AceK) = -1;
S(a.x.Acoa2act,a.f.g.acoa2act) = 1;
S(a.x.Acoa2act,a.f.d.Acoa2act) = -1;
S(a.x.Acs,a.f.g.acs) = 1;
S(a.x.Acs,a.f.d.Acs) = -1;
S(a.x.Akg2mal,a.f.g.akg2mal) = 1;
S(a.x.Akg2mal,a.f.d.Akg2mal) = -1;
S(a.x.CAMPdegr,a.f.g.campdegr) = 1;
S(a.x.CAMPdegr,a.f.d.CAMPdegr) = -1;
S(a.x.Cra,a.f.g.cra) = 1;
S(a.x.Cra,a.f.d.Cra) = -1;
S(a.x.CraFBP,a.f.d.CraFBP) = -1;
S(a.x.Crp,a.f.g.crp) = 1;
S(a.x.Crp,a.f.d.Crp) = -1;
S(a.x.CrpcAMP,a.f.d.CrpcAMP) = -1;
S(a.x.Cya,a.f.g.cya) = 1;
S(a.x.Cya,a.f.d.Cya) = -1;
S(a.x.Emp,a.f.g.emp) = 1;
S(a.x.Emp,a.f.d.Emp) = -1;
S(a.x.Eno,a.f.g.eno) = 1;
S(a.x.Eno,a.f.d.Eno) = -1;
S(a.x.Fdp,a.f.g.fdp) = 1;
S(a.x.Fdp,a.f.d.Fdp) = -1;
S(a.x.GltA,a.f.g.gltA) = 1;
S(a.x.GltA,a.f.d.GltA) = -1;
S(a.x.Icd,a.f.g.icd) = 1;
S(a.x.Icd,a.f.d.Icd) = -1;
S(a.x.Icd_P,a.f.d.Icd_P) = -1;
S(a.x.IclR,a.f.g.iclr) = 1;
S(a.x.IclR,a.f.d.IclR) = -1;
S(a.x.Mdh,a.f.g.mdh) = 1;
S(a.x.Mdh,a.f.d.Mdh) = -1;
S(a.x.MaeAB,a.f.g.maeAB) = 1;
S(a.x.MaeAB,a.f.d.MaeAB) = -1;
S(a.x.PckA,a.f.g.pckA) = 1;
S(a.x.PckA,a.f.d.PckA) = -1;
S(a.x.Pdh,a.f.g.pdh) = 1;
S(a.x.Pdh,a.f.d.Pdh) = -1;
S(a.x.PdhR,a.f.g.pdhr) = 1;
S(a.x.PdhR,a.f.d.PdhR) = -1;
S(a.x.PdhRPYR,a.f.d.PdhRPYR) = -1;
S(a.x.PfkA,a.f.g.pfkA) = 1;
S(a.x.PfkA,a.f.d.PfkA) = -1;
S(a.x.Ppc,a.f.g.ppc) = 1;
S(a.x.Ppc,a.f.d.Ppc) = -1;
S(a.x.PpsA,a.f.g.ppsA) = 1;
S(a.x.PpsA,a.f.d.PpsA) = -1;
S(a.x.PykF,a.f.g.pykF) = 1;
S(a.x.PykF,a.f.d.PykF) = -1;
S(a.x.Icd,a.f.e.AceK_Ki) = -1;
S(a.x.Icd,a.f.e.AceK_Ph) = 1;
S(a.x.Icd_P,a.f.e.AceK_Ki) = 1;
S(a.x.Icd_P,a.f.e.AceK_Ph) = -1;
S(a.x.EIIA,a.f.pts.r1) = -1;
S(a.x.EIIA,a.f.pts.r4) = 1;
S(a.x.EIIA,a.f.g.EIIA) = 1;
S(a.x.EIIA,a.f.d.EIIA) = -1;
S(a.x.EIIA_P,a.f.pts.r1) = 1;
S(a.x.EIIA_P,a.f.pts.r4) = -1;
S(a.x.EIIA_P,a.f.d.EIIA_P) = -1;
S(a.x.EIICB,a.f.g.EIICB) = 1;
S(a.x.EIICB,a.f.d.EIICB) = -1;
S(a.x.CraFBP,a.f.tf.Cra) = 1;
S(a.x.Cra,a.f.tf.Cra) = -1;
S(a.x.CrpcAMP,a.f.tf.Crp) = 1;
S(a.x.Crp,a.f.tf.Crp) = -1;
S(a.x.PdhRPYR,a.f.tf.PdhR) = 1;
S(a.x.PdhR,a.f.tf.PdhR) = -1;
S(a.x.ACoA,a.f.e.Acs) = 1;
S(a.x.ACoA,a.f.e.Pdh) = 1;
S(a.x.ACoA,a.f.e.Acoa2act) = -1;
S(a.x.ACoA,a.f.e.GltA) = -1;
S(a.x.ACoA,a.f.e.AceB) = -1;
S(a.x.ACoA,a.f.d.ACoA) = -1;
S(a.x.ACoA,a.f.bm.ACoA) = -1;
S(a.x.AKG,a.f.e.Icd) = 1;
S(a.x.AKG,a.f.e.AceA) = 1;
S(a.x.AKG,a.f.e.Akg2mal) = -1;
S(a.x.AKG,a.f.d.AKG) = -1;
S(a.x.AKG,a.f.bm.AKG) = -1;
S(a.x.cAMP,a.f.e.Cya) = 1;
S(a.x.cAMP,a.f.e.CAMPdegr) = -1;
S(a.x.cAMP,a.f.d.cAMP) = -1;
S(a.x.FBP,a.f.e.PfkA) = 1;
S(a.x.FBP,a.f.e.Emp) = -0.5;
S(a.x.FBP,a.f.e.Fdp) = -1;
S(a.x.FBP,a.f.d.FBP) = -1;
S(a.x.G6P,a.f.pts.r4) = 1;
S(a.x.G6P,a.f.e.Fdp) = 1;
S(a.x.G6P,a.f.e.PfkA) = -1;
S(a.x.G6P,a.f.d.G6P) = -1;
S(a.x.G6P,a.f.bm.G6P) = -1;
S(a.x.GLX,a.f.e.AceA) = 1;
S(a.x.GLX,a.f.e.AceB) = -1;
S(a.x.GLX,a.f.d.GLX) = -1;
S(a.x.ICT,a.f.e.GltA) = 1;
S(a.x.ICT,a.f.e.AceA) = -1;
S(a.x.ICT,a.f.e.Icd) = -1;
S(a.x.ICT,a.f.d.ICT) = -1;
S(a.x.MAL,a.f.e.AceB) = 1;
S(a.x.MAL,a.f.e.Akg2mal) = 1;
S(a.x.MAL,a.f.e.MaeAB) = -1;
S(a.x.MAL,a.f.e.Mdh) = -1;
S(a.x.MAL,a.f.d.MAL) = -1;
S(a.x.OAA,a.f.e.PckA) = -1;
S(a.x.OAA,a.f.e.GltA) = -1;
S(a.x.OAA,a.f.e.Ppc) = 1;
S(a.x.OAA,a.f.e.Mdh) = 1;
S(a.x.OAA,a.f.d.OAA) = -1;
S(a.x.OAA,a.f.bm.OAA) = -1;
S(a.x.PEP,a.f.e.PckA) = 1;
S(a.x.PEP,a.f.e.PpsA) = 1;
S(a.x.PEP,a.f.e.PykF) = -1;
S(a.x.PEP,a.f.e.Ppc) = -1;
S(a.x.PEP,a.f.e.Eno) = 1;
S(a.x.PEP,a.f.pts.r1) = -1;
S(a.x.PEP,a.f.d.PEP) = -1;
S(a.x.PEP,a.f.bm.PEP) = -1;
S(a.x.PG3,a.f.e.Emp) = 1;
S(a.x.PG3,a.f.e.Eno) = -1;
S(a.x.PG3,a.f.d.PG3) = -1;
S(a.x.PG3,a.f.bm.PG3) = -1;
S(a.x.PYR,a.f.e.MaeAB) = 1;
S(a.x.PYR,a.f.e.PykF) = 1;
S(a.x.PYR,a.f.e.PpsA) = -1;
S(a.x.PYR,a.f.e.Pdh) = -1;
S(a.x.PYR,a.f.pts.r1) = 1;
S(a.x.PYR,a.f.d.PYR) = -1;
S(a.x.PYR,a.f.bm.PYR) = -1;

%% Definition of the parameter values 'p'

p(a.p.env.M_ACT) = 60.05;
p(a.p.env.M_GLC) = 180.156;
p(a.p.env.uc) = 9.5E-07;
p(a.p.e.AceA.kcat) = 614;
p(a.p.e.AceA.n) = 4;
p(a.p.e.AceA.L) = 5.01E+04;
p(a.p.e.AceA.Kict) = 0.022;
p(a.p.e.AceA.Kpep) = 0.055;
p(a.p.e.AceA.Kpg3) = 0.72;
p(a.p.e.AceA.Kakg) = 0.827;
p(a.p.e.AceB.kcat) = 47.8;
p(a.p.e.AceB.Kglx) = 0.95;%0.504;
p(a.p.e.AceB.Kacoa) = 0.755;
p(a.p.e.AceB.Kglxacoa) = 0.719;
p(a.p.e.AceK.kcat_ki) = 3.4E+12;
p(a.p.e.AceK.kcat_ph) = 1.7E+09;
p(a.p.e.AceK.n) = 2;
p(a.p.e.AceK.L) = 1.0E+08;
p(a.p.e.AceK.Kicd) = 0.043;
p(a.p.e.AceK.Kicd_P) = 0.643;
p(a.p.e.AceK.Kpep) = 0.539;
p(a.p.e.AceK.Kpyr) = 0.038;
p(a.p.e.AceK.Koaa) = 0.173;
p(a.p.e.AceK.Kglx) = 0.866;
p(a.p.e.AceK.Kakg) = 0.82;
p(a.p.e.AceK.Kpg3) = 1.57;
p(a.p.e.AceK.Kict) = 0.137;
p(a.p.e.Acoa2act.kcat) = 3079;
p(a.p.e.Acoa2act.n) = 2;
p(a.p.e.Acoa2act.L) = 6.39E+05;
p(a.p.e.Acoa2act.Kacoa) = 0.022;
p(a.p.e.Acoa2act.Kpyr) = 0.022;
p(a.p.e.Acs.kcat) = 340/0.000036201*0.001096222;
p(a.p.e.Acs.Kact) = 1.0E-03;
p(a.p.e.Akg2mal.kcat) = 1530;
p(a.p.e.Akg2mal.Kakg) = 0.548;
p(a.p.e.CAMPdegr.kcat) = 1.00E+03;
p(a.p.e.CAMPdegr.KcAMP) = 0.1;
p(a.p.e.Cya.kcat) = 993;
p(a.p.e.Cya.KEIIA) = 1.7E-03;
p(a.p.e.Emp.kcat.f) = 1000/0.011389032*0.011515593;
p(a.p.e.Emp.kcat.r) = 848/0.011389032*0.011515593;
p(a.p.e.Emp.Kfbp) = 5.92;
p(a.p.e.Emp.Kpg3) = 16.6;
p(a.p.e.Eno.kcatf) = 695/0.011389032*0.011552813;
p(a.p.e.Eno.kcatr) = 522/0.011389032*0.011552813;
p(a.p.e.Eno.Kpg3) = 4.76;
p(a.p.e.Eno.Kpep) = 1.11;
p(a.p.e.Fdp.kcat) = 192/0.000074810*0.000157492;
p(a.p.e.Fdp.n) = 4;
p(a.p.e.Fdp.L) = 4.0E+06;
p(a.p.e.Fdp.Kfbp) = 3.0E-03;
p(a.p.e.Fdp.Kpep) = 0.3;
p(a.p.e.GltA.kcat) = 1614/0.000292771*0.001029612;
p(a.p.e.GltA.Koaa) = 0.029;
p(a.p.e.GltA.Kacoa) = 0.212;
p(a.p.e.GltA.Koaaacoa) = 0.029;
p(a.p.e.GltA.Kakg) = 0.63;
p(a.p.e.Icd.kcat) = 695;
p(a.p.e.Icd.n) = 2;
p(a.p.e.Icd.L) = 127;
p(a.p.e.Icd.Kict) = 1.6E-04;
p(a.p.e.Icd.Kpep) = 0.334;
p(a.p.e.Mdh.kcat) = 773/0.000491491*0.00345727;
p(a.p.e.Mdh.n) = 1.7;
p(a.p.e.Mdh.Kmal) = 10.1;
p(a.p.e.MaeAB.kcat) = 1879;
p(a.p.e.MaeAB.n) = 1.33;
p(a.p.e.MaeAB.L) = 1.04E+05;
p(a.p.e.MaeAB.Kmal) = 6.24E-03;
p(a.p.e.MaeAB.Kacoa) = 3.64;
p(a.p.e.MaeAB.Kcamp) = 6.54;
p(a.p.e.PckA.kcat) = 55.5/0.000336947*0.002290892;
p(a.p.e.PckA.Koaa) = 0.184;
p(a.p.e.PckA.Kpep) = 1000;
p(a.p.e.Pdh.kcat) = 1179/0.001*0.004647401;
p(a.p.e.Pdh.n) = 2.65;
p(a.p.e.Pdh.L) = 3.4;
p(a.p.e.Pdh.Kpyr) = 0.128;
p(a.p.e.Pdh.KpyrI) = 0.231;
p(a.p.e.Pdh.Kglx) = 0.218;
p(a.p.e.PfkA.kcat) = 9.08E+05/0.000242131*0.000143816;
p(a.p.e.PfkA.n) = 4;
p(a.p.e.PfkA.L) = 9.5E+07;
p(a.p.e.PfkA.Kg6p) = 0.022;
p(a.p.e.PfkA.Kpep) = 0.138;
p(a.p.e.Ppc.kcat) = 5635/0.000377962*0.000999714;
p(a.p.e.Ppc.n) = 3;
p(a.p.e.Ppc.L) = 5.2E+06;
p(a.p.e.Ppc.Kpep) = 0.048;
p(a.p.e.Ppc.Kfbp) = 0.408;
p(a.p.e.PpsA.kcat) = 1.32;
p(a.p.e.PpsA.n) = 2;
p(a.p.e.PpsA.L) = 1.0E-79;
p(a.p.e.PpsA.Kpyr) = 1.77E-03;
p(a.p.e.PpsA.Kpep) = 1.0E-03;
p(a.p.e.PykF.kcat) = 5749/0.002501893*0.005977168;
p(a.p.e.PykF.n) = 4;
p(a.p.e.PykF.L) = 1.0E+05;
p(a.p.e.PykF.Kpep) = 5;
p(a.p.e.PykF.Kfbp) = 0.413;
p(a.p.pts.k1) = 116;
p(a.p.pts.km1) = 46.3;
p(a.p.pts.k4) = 2520;
p(a.p.pts.KEIIA) = 8.5E-03;
p(a.p.pts.Kglc) = 1.2E-03;
p(a.p.tf.Cra.scale) = 100;
p(a.p.tf.Cra.kfbp) = 1.36;
p(a.p.tf.Cra.n) = 2;
p(a.p.tf.Crp.scale) = 1.0E+08;
p(a.p.tf.Crp.kcamp) = 0.895;
p(a.p.tf.Crp.n) = 1;
p(a.p.tf.PdhR.scale) = 100;
p(a.p.tf.PdhR.kpyr) = 0.164;
p(a.p.tf.PdhR.n) = 1;
p(a.p.g.aceBAK.vcra_unbound) = 1.9E-09;
p(a.p.g.aceBAK.vcra_bound) = 2.0E-06;
p(a.p.g.aceBAK.Kcra) = 3.65E-03;
p(a.p.g.aceBAK.aceBfactor) = 0.3;
p(a.p.g.aceBAK.aceKfactor) = 0.03;
p(a.p.g.aceBAK.KDNA) = 2.19;
p(a.p.g.aceBAK.KP) = 0.897;
p(a.p.g.aceBAK.KPprime) = 3.01E-03;
p(a.p.g.aceBAK.KG) = 4.88E-03;
p(a.p.g.aceBAK.L) = 923;
p(a.p.g.aceBAK.kcat_iclr) = 9.3E-04;
p(a.p.g.aceBAK.DNA) = 1;
p(a.p.g.aceBAK.vcrp_bound) = 2.3E-10;
p(a.p.g.aceBAK.vcrp_unbound) = 2.0E-08;
p(a.p.g.aceBAK.Kcrp) = 0.341;
p(a.p.g.acs.vcrp_unbound) = 0;
p(a.p.g.acs.vcrp_bound) = 1.2E-06*0.000036201/0.001096222;
p(a.p.g.acs.n) = 2.31;
p(a.p.g.acs.Kcrp) = 4.7E-03;
p(a.p.g.akg2mal.vcrp_unbound) = 0;
p(a.p.g.akg2mal.vcrp_bound) = 1.4E-06;
p(a.p.g.akg2mal.Kcrp) = 0.091;
p(a.p.g.akg2mal.n) = 0.74;
p(a.p.g.emp.vcra_unbound) = 6.2E-07*0.011389032/0.011515593;
p(a.p.g.emp.vcra_bound) = 0;
p(a.p.g.emp.Kcra) = 0.09;
p(a.p.g.emp.vcrp_unbound) = 0;
p(a.p.g.emp.vcrp_bound) = 4.7E-07;
p(a.p.g.emp.Kcrp) = 0.012;
p(a.p.g.eno.vcra_unbound) = 6.8E-07*0.011389032/0.011552813;
p(a.p.g.eno.vcra_bound) = 0;
p(a.p.g.eno.Kcra) = 0.016;
p(a.p.g.fdp.vcra_unbound) = 0;
p(a.p.g.fdp.vcra_bound) = 4.5E-08*0.000074810/0.000157492;
p(a.p.g.fdp.Kcra) = 1.18E-03;
p(a.p.g.gltA.vcrp_unbound) = 0;
p(a.p.g.gltA.vcrp_bound) = 2.3E-06*0.000292771/0.001029612;
p(a.p.g.gltA.Kcrp) = 0.04;
p(a.p.g.gltA.n) = 1.07;
p(a.p.g.icd.vcra_unbound) = 1.1E-07;
p(a.p.g.icd.vcra_bound) = 8.5E-07;
p(a.p.g.icd.Kcra) = 1.17E-03;
p(a.p.g.mdh.vcrp_unbound) = 0;
p(a.p.g.mdh.vcrp_bound) = 9.1E-06*0.000491491/0.00345727;
p(a.p.g.mdh.Kcrp) = 0.06;
p(a.p.g.pckA.vcra_unbound) = 0;
p(a.p.g.pckA.vcra_bound) = 2.5E-06*0.000336947/0.002290892;
p(a.p.g.pckA.Kcra) = 5.35E-03;
p(a.p.g.pdh.vpdhr_unbound) = 3.6E-07*0.001/0.004647401;
p(a.p.g.pdh.vpdhr_bound) = 1.3E-09*0.001/0.004647401;
p(a.p.g.pdh.Kpdhr) = 3.4E-03;
p(a.p.g.pfkA.vcra_unbound) = 8.2E-07*0.000242131/0.000143816;
p(a.p.g.pfkA.vcra_bound) = 6.6E-09*0.000242131/0.000143816;
p(a.p.g.pfkA.Kcra) = 6.3E-07;
p(a.p.g.ppsA.vcra_unbound) = 0;
p(a.p.g.ppsA.vcra_bound) = 3.3E-06;
p(a.p.g.ppsA.Kcra) = 0.017;
p(a.p.g.pykF.vcra_unbound) = 3.9E-07*0.002501893/0.005977168;
p(a.p.g.pykF.vcra_bound) = 2.1E-09*0.002501893/0.005977168;
p(a.p.g.pykF.Kcra) = 2.3E-03;
p(a.p.d.k_degr) = 2.8E-05;
p(a.p.bm.k_expr) = 2.0E+04;
p(a.p.bm.muACT) = 5.6E-05;
p(a.p.bm.muGLC) = 1.8E-04;
p(a.p.bm.GLC.ACoA) = 1.88;
p(a.p.bm.GLC.AKG) = 0.978;
p(a.p.bm.GLC.G6P) = 0.154;
p(a.p.bm.GLC.OAA) = 6.4;
p(a.p.bm.GLC.PEP) = 0.423;
p(a.p.bm.GLC.PG3) = 0.049;
p(a.p.bm.GLC.PYR) = 0.553;
p(a.p.bm.ACT.ACoA) = 0.108;
p(a.p.bm.ACT.AKG) = 0.056;
p(a.p.bm.ACT.G6P) = 0.076;
p(a.p.bm.ACT.OAA) = 1.43;
p(a.p.bm.ACT.PEP) = 0.047;
p(a.p.bm.ACT.PG3) = 0.066;
p(a.p.bm.ACT.PYR) = 5.185;

%% Write steady state values into auxiliary variables
% The two auxiliary variables ssACT and ssGLC contain the
% steady states on acetate and glucose, respectively.
% Remember that these values, which have been obtained with
% simulations on acetate and glucose as sole carbon sources,
% depend on the chosen parameter set.
% The sole purpose of these two auxiliary variables is to allow
% for the initial conditions to be set to either steady state.
% An initial condition equal to a steady state is interesting
% when a transition to the other steady state is investigated.

ssACT(a.x.ACoA) = 2.045970006;
ssACT(a.x.AKG) = 1.11161563;
ssACT(a.x.cAMP) = 4.021170806;
ssACT(a.x.FBP) = 0.272276143;
ssACT(a.x.G6P) = 1.145931991;
ssACT(a.x.GLX) = 1.322800153;
ssACT(a.x.ICT) = 1.480239304;
ssACT(a.x.MAL) = 6.396156035;
ssACT(a.x.OAA) = 0.064573174;
ssACT(a.x.PEP) = 0.557008135;
ssACT(a.x.PG3) = 1.293579693;
ssACT(a.x.PYR) = 0.037723095;
ssACT(a.x.AceA) = 0.101787757;
ssACT(a.x.AceB) = 0.030536327;
ssACT(a.x.AceK) = 0.003053633;
ssACT(a.x.Acoa2act) = 0.001;
ssACT(a.x.Acs) = 0.010144567*0.000036201/0.001096222;
ssACT(a.x.Akg2mal) = 0.002192398;
ssACT(a.x.CAMPdegr) = 0.001;
ssACT(a.x.Cya) = 0.001;
ssACT(a.x.Emp) = 0.009751833*0.011389032/0.011515593;
ssACT(a.x.Eno) = 0.006304314*0.011389032/0.011552813;
ssACT(a.x.Fdp) = 0.000513512*0.000074810/0.000157492;
ssACT(a.x.GltA) = 0.003539257*0.000292771/0.001029612;
ssACT(a.x.Icd) = 0.002404566;
ssACT(a.x.Icd_P) = 0.007516755;
ssACT(a.x.Mdh) = 0.010969029*0.000491491/0.00345727;
ssACT(a.x.MaeAB) = 0.003399346;
ssACT(a.x.PckA) = 0.018902966*0.000336947/0.002290892;
ssACT(a.x.Pdh) = 0.001760705*0.001/0.004647401;
ssACT(a.x.PfkA) = 8.89703E-05*0.000242131/0.000143816;
ssACT(a.x.Ppc) = 0.000279893*0.000377962/0.000999714;
ssACT(a.x.PpsA) = 0.012844496;
ssACT(a.x.PykF) = 0.001305745*0.002501893/0.005977168;
ssACT(a.x.Cra) = 0.007009039;
ssACT(a.x.CraFBP) = 0.000280931;
ssACT(a.x.Crp) = 0.001327161;
ssACT(a.x.CrpcAMP) = 0.005962839;
ssACT(a.x.IclR) = 0.00729;
ssACT(a.x.PdhR) = 0.005926738;
ssACT(a.x.PdhRPYR) = 0.001363262;
ssACT(a.x.EIIA) = 0.002631995;
ssACT(a.x.EIIA_P) = 0.097368005;
ssACT(a.x.EIICB) = 0.003;

ssGLC(a.x.ACoA) = 0.351972298;
ssGLC(a.x.AKG) = 0.191190619;
ssGLC(a.x.cAMP) = 0.202804098;
ssGLC(a.x.FBP) = 6.57504207;
ssGLC(a.x.G6P) = 1.908140784;
ssGLC(a.x.GLX) = 5.70593E-09;
ssGLC(a.x.ICT) = 0.001408116;
ssGLC(a.x.MAL) = 3.278779135;
ssGLC(a.x.OAA) = 0.050535354;
ssGLC(a.x.PEP) = 0.210455879;
ssGLC(a.x.PG3) = 5.720977255;
ssGLC(a.x.PYR) = 0.863278018;
ssGLC(a.x.AceA) = 0.00472323;
ssGLC(a.x.AceB) = 0.001416969;
ssGLC(a.x.AceK) = 0.000141697;
ssGLC(a.x.Acoa2act) = 0.001;
ssGLC(a.x.Acs) = 0.000036201;%old: 0.001096222;
ssGLC(a.x.Akg2mal) = 0.001026848;
ssGLC(a.x.CAMPdegr) = 0.001;
ssGLC(a.x.Cya) = 0.001;
ssGLC(a.x.Emp) = 0.011389032;%old: 0.011515593;
ssGLC(a.x.Eno) = 0.011389032;%old: 0.011552813;
ssGLC(a.x.Fdp) = 0.000074810;%old: 0.000157492;
ssGLC(a.x.GltA) = 0.000292771;%old: 0.001029612;
ssGLC(a.x.Icd) = 0.004290789;
ssGLC(a.x.Icd_P) = 0.000220477;
ssGLC(a.x.Mdh) = 0.000491491;% old: 0.00345727;
ssGLC(a.x.MaeAB) = 0.000999714;
ssGLC(a.x.PckA) = 0.000336947;%old: 0.002290892;
ssGLC(a.x.Pdh) = 0.001;% old: 0.004647401;
ssGLC(a.x.PfkA) = 0.000242131;%old: 0.000143816;
ssGLC(a.x.Ppc) = 0.000377962;%old: 0.000999714;
ssGLC(a.x.PpsA) = 0.000987493;
ssGLC(a.x.PykF) = 0.002501893;%old: 0.005977168;
ssGLC(a.x.Cra) = 0.000299098;
ssGLC(a.x.CraFBP) = 0.006990902;
ssGLC(a.x.Crp) = 0.005943273;
ssGLC(a.x.CrpcAMP) = 0.001346727;
ssGLC(a.x.IclR) = 0.00729;
ssGLC(a.x.PdhR) = 0.001163813;
ssGLC(a.x.PdhRPYR) = 0.006126187;
ssGLC(a.x.EIIA) = 0.09647707;
ssGLC(a.x.EIIA_P) = 0.00352292;
ssGLC(a.x.EIICB) = 0.003;


%% Set simulation time and initial conditions
% These values are set according to the scenario chosen at the very
% beginning of this file.
% If the scenario 7 'custom' is chosen, the desired values can be entered
% in that section.

% Scenario 1: The population is adapted to glucose,
% the carbon source is glucose

if scenario==1,
    t0 = 0;
    tf = 10 * 3600;
    ic = ssGLC;
    ic(a.x.ACT) = 0;
    ic(a.x.GLC) = 5;
    ic(a.x.OD) = 1.0E-03;
end

% Scenario 2: The population is adapted to acetate,
% the carbon source is acetate

if scenario==2,
    t0 = 0;
    tf = 10 * 3600;
    ic = ssACT;
    ic(a.x.ACT) = 5;
    ic(a.x.GLC) = 0;
    ic(a.x.OD) = 1.0E-03;
end   

% Scenario 3: The population is adapted to glucose,
% the carbon source is acetate

if scenario==3,
    t0 = 0;
    tf = 30 * 3600;
    ic = ssGLC;
    ic(a.x.ACT) = 5;
    ic(a.x.GLC) = 0;
    ic(a.x.OD) = 1.0E-04;
end 

% Scenario 4: The population is adapted to acetate,
% the carbon source is glucose

if scenario==4,
    t0 = 0;
    tf = 10 * 3600;
    ic = ssACT;
    ic(a.x.ACT) = 0;
    ic(a.x.GLC) = 5;
    ic(a.x.OD) = 1.0E-04;
end

% Scenario 5: Diauxic shift scenario as presented in the paper

if scenario==5,
    t0 = 0;
    tf = 10 * 3600;
    ic = ssGLC;
    ic(a.x.ACT) = 0;
    ic(a.x.GLC) = 4.8;
    ic(a.x.OD) = 0.03;
end

% Scenario 6: Extended diauxic shift scenario as presented in the paper
% As this scenario is more complex, this special case is further dealt with
% below.

if scenario==6,
    t0 = 0;
    tf = 8.15 * 3600; %6.7
    ic = ssGLC;
    ic(a.x.ACT) = 0;
    ic(a.x.GLC) = 4.8;
    ic(a.x.OD) = 0.03;
end

% Scenario 7: Custom simulation
% To run a custom simulation, set scenario=7 at the very beginning of this
% file, and enter the desired simulation time and initial
% conditions here. To change parameter values, search for the string
% 'Definition of the parameter values' and change the values there.

if scenario==7,
    
    % Enter the initial time here, in seconds
    t0 = 0;
    
    % Enter the final time here, in seconds
    tf = 10 * 3600;
    
    % Enter the initial conditions of the intracellular states here
    % To enter the initial conditions, you can use
    % ic = ssACT; % to set the I.C.'s to the steady state on acetate, or
     ic = ssGLC; % to set the I.C.'s to the steady state on glucose
    % Alternatively or in combination with this, you can define the
    % initial conditions of individual states through, for instance, 
    % ic(a.x.ACoA) = 0.349694268;
    % ic(a.x.AKG) = 0.195376075;
    
    % Enter the initial acetate concentration here, in g/l
    ic(a.x.ACT) = 0; % g/l
    
    % Enter the initial glucose concentration here, in g/l
    ic(a.x.GLC) = 4.8; % g/l
    
    % Enter the initial optical density here, in OD_600
    ic(a.x.OD)=0.03;
    
end


%% Run the simulation

options=odeset('RelTol',1e-5,'AbsTol',1e-20);

called_by_ode15s = true; % tell ecoliOK_core to return state derivatives dx
[t x] = ode15s(@ecoliOK_core,[t0 tf],ic,options); % run simulation
    
if scenario==6, % Append the second simulation run of scenario 6
    
    % Inoculate a 50ml shake flask containing 5 g/l acetate using cells
    % from the previous culture
    
    new_ic = x(end,:);
    new_ic(a.x.ACT) = 5;
    new_ic(a.x.GLC) = 0;
    new_ic(a.x.OD) = 0.03; % 0.02, 21.6
    new_tf = (t(end)/3600 + 19.7) * 3600;

    % Run the second simulation
    
    [t2 x2] = ode15s(@ecoliOK_core,[t(end) new_tf],new_ic,options);
    
    % Inoculate a 50ml shake flask containing 5 g/l acetate and glucose
    % using cells from the previous culture
    
    new_ic = x2(end,:);
    new_ic(a.x.ACT) = 3;
    new_ic(a.x.GLC) = 3;
    new_ic(a.x.OD) = 5E-4; 
    new_tf = (t2(end)/3600 + 16.45) * 3600; % 3 3 5E-4 16.45 ok

    % Run the third simulation
    
    [t3 x3] = ode15s(@ecoliOK_core,[t2(end) new_tf],new_ic,options);
    
    % Fuse the results of the three simulation runs
    
    x=[x; x2; x3];
    t=[t; t2; t3];
end

%% Write the simulation results into ecoliOK's output structure 'sim'

sim.t = t;
sim.x = x;

called_by_ode15s = false; % tell ecoliOK_core to return reaction rates f
for indx = 1:size(x,1), % reconstruct the reaction rates f from the states x
    sim.f(indx,:) = ecoliOK_core(0,x(indx,:));
end

sim.p = p;
sim.ic = ic;
sim.S = S;
sim.a = a;
sim.scenario = scenario;

% This is the end of the main function code.
% The following code concerns the nested function ecoliOK_core.


%% ODE model
% The ODE model is contained in the nested function ecoliOK_core,
% and is called by the main function's ode15s solver.
% ecoliOK_core receives the state variables x (and the time t) from
% the ode15s solver, and returns the change of x over time, dx,
% to the solver.
% ecoliOK_core reads the variables S (stoichiometric matrix),
% a (aliases), p (parameters), and the boolean flag called_by_ode15s,
% which were defined in the main function.
    
    function vargout = ecoliOK_core(t,x)
  
        %% Growth rate, biomass reactions, and dependent enzyme levels
        % As a simplification, in this model, the growth rate,
        % the seven biomass reaction rates, and the actual steady state
        % levels of the two enzymes MaeAB and Ppc depend on the carbon
        % sources present in the environment.
        
        % Interpolate between pure glucose and pure acetate
        % growth conditions
        alphaGLC = x(a.x.GLC)/(x(a.x.GLC)+p(a.p.pts.Kglc));
        alphaACT = x(a.x.ACT)/(x(a.x.ACT)+p(a.p.e.Acs.Kact))*(1-x(a.x.GLC)/(x(a.x.GLC)+p(a.p.pts.Kglc)));
  
        % Calculate the growth rate 'mu'
        mu = alphaGLC*p(a.p.bm.muGLC) + alphaACT*p(a.p.bm.muACT);
        
        % Calculate the first order rate constants of the seven biomass reactions
        k_bm_ACoA = alphaGLC*p(a.p.bm.GLC.ACoA) + alphaACT*p(a.p.bm.ACT.ACoA);
        k_bm_AKG = alphaGLC*p(a.p.bm.GLC.AKG) + alphaACT*p(a.p.bm.ACT.AKG);       
        k_bm_G6P = alphaGLC*p(a.p.bm.GLC.G6P) + alphaACT*p(a.p.bm.ACT.G6P);
        k_bm_OAA = alphaGLC*p(a.p.bm.GLC.OAA) + alphaACT*p(a.p.bm.ACT.OAA);
        k_bm_PEP = alphaGLC*p(a.p.bm.GLC.PEP) + alphaACT*p(a.p.bm.ACT.PEP);
        k_bm_PG3 = alphaGLC*p(a.p.bm.GLC.PG3) + alphaACT*p(a.p.bm.ACT.PG3);
        k_bm_PYR = alphaGLC*p(a.p.bm.GLC.PYR) + alphaACT*p(a.p.bm.ACT.PYR);
        
        % Calculate biomass reaction rates with 1st order kinetics
        f([
        a.f.bm.ACoA
        a.f.bm.AKG
        a.f.bm.G6P
        a.f.bm.OAA
        a.f.bm.PEP
        a.f.bm.PG3
        a.f.bm.PYR
        ]) = [
        k_bm_ACoA
        k_bm_AKG
        k_bm_G6P
        k_bm_OAA
        k_bm_PEP
        k_bm_PG3
        k_bm_PYR
        ] .* [
        x(a.x.ACoA)
        x(a.x.AKG)
        x(a.x.G6P)
        x(a.x.OAA)
        x(a.x.PEP)
        x(a.x.PG3)
        x(a.x.PYR) ];
    
        % Calculate the actual steady state levels of MaeAB and Ppc
        SS_MaeAB = alphaGLC*ssGLC(a.x.MaeAB) + alphaACT*ssACT(a.x.MaeAB);
        SS_Ppc = alphaGLC*ssGLC(a.x.Ppc) + alphaACT*ssACT(a.x.Ppc);
     
        
        %% Protein phosphorylation rates
        
        % PTS phosphorylation kinetics
        f(a.f.pts.r1) = p(a.p.pts.k1)*x(a.x.PEP)*x(a.x.EIIA)-p(a.p.pts.km1)*x(a.x.PYR)*x(a.x.EIIA_P);
        f(a.f.pts.r4) = p(a.p.pts.k4)*x(a.x.EIICB)*x(a.x.EIIA_P)*x(a.x.GLC)/((p(a.p.pts.KEIIA)+x(a.x.EIIA_P))*(p(a.p.pts.Kglc)+x(a.x.GLC)));
        
        % AceK_ki kinetics: MWC, substrate: Icd, inhibitors: GLX, ICT, OAA,
        % PEP, PG3, PYR, AKG
        f(a.f.e.AceK_Ki) = x(a.x.AceK)*p(a.p.e.AceK.kcat_ki)*x(a.x.Icd)/p(a.p.e.AceK.Kicd)*(1+x(a.x.Icd)/p(a.p.e.AceK.Kicd))^(p(a.p.e.AceK.n)-1)/((1+x(a.x.Icd)/p(a.p.e.AceK.Kicd))^p(a.p.e.AceK.n)+p(a.p.e.AceK.L)*(1+x(a.x.ICT)/p(a.p.e.AceK.Kict)+x(a.x.GLX)/p(a.p.e.AceK.Kglx)+x(a.x.OAA)/p(a.p.e.AceK.Koaa)+x(a.x.AKG)/p(a.p.e.AceK.Kakg)+x(a.x.PEP)/p(a.p.e.AceK.Kpep)+x(a.x.PG3)/p(a.p.e.AceK.Kpg3)+x(a.x.PYR)/p(a.p.e.AceK.Kpyr))^p(a.p.e.AceK.n));
        
        % AceK_ph kinetics: MWC, substrate: Icd_P, activators: OAA, PEP,
        % PG3, PYR, AKG
        f(a.f.e.AceK_Ph) = x(a.x.AceK)*p(a.p.e.AceK.kcat_ph)*x(a.x.Icd_P)/p(a.p.e.AceK.Kicd_P)*(1+x(a.x.Icd_P)/p(a.p.e.AceK.Kicd_P))^(p(a.p.e.AceK.n)-1)/((1+x(a.x.Icd_P)/p(a.p.e.AceK.Kicd_P))^p(a.p.e.AceK.n)+p(a.p.e.AceK.L)/(1+x(a.x.OAA)/p(a.p.e.AceK.Koaa)+x(a.x.AKG)/p(a.p.e.AceK.Kakg)+x(a.x.PEP)/p(a.p.e.AceK.Kpep)+x(a.x.PG3)/p(a.p.e.AceK.Kpg3)+x(a.x.PYR)/p(a.p.e.AceK.Kpyr))^p(a.p.e.AceK.n)); 
 
        
        %% Metabolite- transcription factor binding rates
        % The IclR-GLX-PYR binding state is incorporated into the gene
        % expression kinetics of the aceBAK operon
        
        % Cra-FBP binding kinetics: Hill
        f(a.f.tf.Cra) = p(a.p.tf.Cra.scale)*((x(a.x.Cra)+x(a.x.CraFBP))*x(a.x.FBP)^p(a.p.tf.Cra.n)/(x(a.x.FBP)^p(a.p.tf.Cra.n)+p(a.p.tf.Cra.kfbp)^p(a.p.tf.Cra.n))-x(a.x.CraFBP));
        
        % Crp-cAMP binding kinetics: Hill
        f(a.f.tf.Crp) = p(a.p.tf.Crp.scale)*((x(a.x.Crp)+x(a.x.CrpcAMP))*x(a.x.cAMP)^p(a.p.tf.Crp.n)/(x(a.x.cAMP)^p(a.p.tf.Crp.n)+p(a.p.tf.Crp.kcamp)^p(a.p.tf.Crp.n))-x(a.x.CrpcAMP));
        
        % PdhR-PYR binding kinetics: Hill
        f(a.f.tf.PdhR) = p(a.p.tf.PdhR.scale)*((x(a.x.PdhR)+x(a.x.PdhRPYR))*x(a.x.PYR)^p(a.p.tf.PdhR.n)/(x(a.x.PYR)^p(a.p.tf.PdhR.n)+p(a.p.tf.PdhR.kpyr)^p(a.p.tf.PdhR.n))-x(a.x.PdhRPYR));

        
        %% Metabolic reaction rates
        
        % AceA kinetics: MWC, substrate: ICT, inhibitors: PG3, PEP, AKG      
        f(a.f.e.AceA) = x(a.x.AceA)*p(a.p.e.AceA.kcat)*x(a.x.ICT)/p(a.p.e.AceA.Kict)*(1+x(a.x.ICT)/p(a.p.e.AceA.Kict))^(p(a.p.e.AceA.n)-1)/((1+x(a.x.ICT)/p(a.p.e.AceA.Kict))^p(a.p.e.AceA.n)+p(a.p.e.AceA.L)*(1+x(a.x.PEP)/p(a.p.e.AceA.Kpep)+x(a.x.PG3)/p(a.p.e.AceA.Kpg3)+x(a.x.AKG)/p(a.p.e.AceA.Kakg))^p(a.p.e.AceA.n));
        
        % AceB kinetics: Two-substrate MM, substrates: GLX, ACoA   
        f(a.f.e.AceB) = x(a.x.AceB)*p(a.p.e.AceB.kcat)*x(a.x.GLX)*x(a.x.ACoA)/(p(a.p.e.AceB.Kglxacoa)*p(a.p.e.AceB.Kacoa)+p(a.p.e.AceB.Kacoa)*x(a.x.GLX)+p(a.p.e.AceB.Kglx)*x(a.x.ACoA)+x(a.x.GLX)*x(a.x.ACoA));
           
        % Acoa2act kinetics: MWC, substrate: ACoA, activator: PYR
        f(a.f.e.Acoa2act) = x(a.x.Acoa2act)*p(a.p.e.Acoa2act.kcat)*x(a.x.ACoA)/p(a.p.e.Acoa2act.Kacoa)*(1+x(a.x.ACoA)/p(a.p.e.Acoa2act.Kacoa))^(p(a.p.e.Acoa2act.n)-1)/((1+x(a.x.ACoA)/p(a.p.e.Acoa2act.Kacoa))^p(a.p.e.Acoa2act.n)+p(a.p.e.Acoa2act.L)/(1+x(a.x.PYR)/p(a.p.e.Acoa2act.Kpyr))^p(a.p.e.Acoa2act.n));
        
        % Acs kinetics: MM, substrate: ACT
        f(a.f.e.Acs) = x(a.x.Acs)*p(a.p.e.Acs.kcat)*x(a.x.ACT)/(x(a.x.ACT)+p(a.p.e.Acs.Kact));
        
        % Akg2mal kinetics: MM, substrate: AKG
        f(a.f.e.Akg2mal) = x(a.x.Akg2mal)*p(a.p.e.Akg2mal.kcat)*x(a.x.AKG)/(x(a.x.AKG)+p(a.p.e.Akg2mal.Kakg));
        
        % CAMPdegr kinetics: MM, substrate: cAMP
        f(a.f.e.CAMPdegr) = p(a.p.e.CAMPdegr.kcat)*x(a.x.CAMPdegr)*x(a.x.cAMP)/(x(a.x.cAMP)+p(a.p.e.CAMPdegr.KcAMP));
        
        % Cya kinetics: MM, substrate: Cya
        f(a.f.e.Cya) = p(a.p.e.Cya.kcat)*x(a.x.Cya)*x(a.x.EIIA_P)/(x(a.x.EIIA_P)+p(a.p.e.Cya.KEIIA));
        
        % Emp kinetics: reversible MM, substrates: FBP, PG3
        f(a.f.e.Emp) = (x(a.x.Emp)*p(a.p.e.Emp.kcat.f)*x(a.x.FBP)/p(a.p.e.Emp.Kfbp)-x(a.x.Emp)*p(a.p.e.Emp.kcat.r)*x(a.x.PG3)/p(a.p.e.Emp.Kpg3))/(1+x(a.x.FBP)/p(a.p.e.Emp.Kfbp)+x(a.x.PG3)/p(a.p.e.Emp.Kpg3));
        
        % Eno kinetics: reversible MM, substrates: PG3, PEP
        f(a.f.e.Eno) = (x(a.x.Eno)*p(a.p.e.Eno.kcatf)*x(a.x.PG3)/p(a.p.e.Eno.Kpg3)-x(a.x.Eno)*p(a.p.e.Eno.kcatr)*x(a.x.PEP)/p(a.p.e.Eno.Kpep))/(1+x(a.x.PG3)/p(a.p.e.Eno.Kpg3)+x(a.x.PEP)/p(a.p.e.Eno.Kpep));
        
        % Fdp kinetics: MWC, substrate: FBP, activator: PEP
        f(a.f.e.Fdp) = x(a.x.Fdp)*p(a.p.e.Fdp.kcat)*x(a.x.FBP)/p(a.p.e.Fdp.Kfbp)*(1+x(a.x.FBP)/p(a.p.e.Fdp.Kfbp))^(p(a.p.e.Fdp.n)-1)/((1+x(a.x.FBP)/p(a.p.e.Fdp.Kfbp))^p(a.p.e.Fdp.n)+p(a.p.e.Fdp.L)/(1+x(a.x.PEP)/p(a.p.e.Fdp.Kpep))^p(a.p.e.Fdp.n));
        
        % GltA kinetics: Two-substrate MM, substrates: OAA, ACoA,
        % competitive inhibitor: AKG
        f(a.f.e.GltA) = x(a.x.GltA)*p(a.p.e.GltA.kcat)*x(a.x.OAA)*x(a.x.ACoA)/((1+x(a.x.AKG)/p(a.p.e.GltA.Kakg))*p(a.p.e.GltA.Koaaacoa)*p(a.p.e.GltA.Kacoa)+p(a.p.e.GltA.Kacoa)*x(a.x.OAA)+(1+x(a.x.AKG)/p(a.p.e.GltA.Kakg))*p(a.p.e.GltA.Koaa)*x(a.x.ACoA)+x(a.x.OAA)*x(a.x.ACoA));
        
        % Icd kinetics: MWC, substrate: ICT, inhibitor: PEP
        f(a.f.e.Icd) = x(a.x.Icd)*p(a.p.e.Icd.kcat)*x(a.x.ICT)/p(a.p.e.Icd.Kict)*(1+x(a.x.ICT)/p(a.p.e.Icd.Kict))^(p(a.p.e.Icd.n)-1)/((1+x(a.x.ICT)/p(a.p.e.Icd.Kict))^p(a.p.e.Icd.n)+p(a.p.e.Icd.L)*(1+x(a.x.PEP)/p(a.p.e.Icd.Kpep))^p(a.p.e.Icd.n));
        
        % MaeAB kinetics: MWC, substrate: MAL, inhibitors: AcoA, cAMP
        f(a.f.e.MaeAB) = x(a.x.MaeAB)*p(a.p.e.MaeAB.kcat)*x(a.x.MAL)/p(a.p.e.MaeAB.Kmal)*(1+x(a.x.MAL)/p(a.p.e.MaeAB.Kmal))^(p(a.p.e.MaeAB.n)-1)/((1+x(a.x.MAL)/p(a.p.e.MaeAB.Kmal))^p(a.p.e.MaeAB.n)+p(a.p.e.MaeAB.L)*(1+x(a.x.ACoA)/p(a.p.e.MaeAB.Kacoa)+x(a.x.cAMP)/p(a.p.e.MaeAB.Kcamp))^p(a.p.e.MaeAB.n));
        
        % Mdh kinetics: Hill, substrate: MAL
        f(a.f.e.Mdh) = x(a.x.Mdh)*p(a.p.e.Mdh.kcat)*x(a.x.MAL)^p(a.p.e.Mdh.n)/(x(a.x.MAL)^p(a.p.e.Mdh.n)+p(a.p.e.Mdh.Kmal)^p(a.p.e.Mdh.n));
        
        % PckA kinetics: MM, substrate: OAA, competitive inhibitor: PEP
        f(a.f.e.PckA) = x(a.x.PckA)*p(a.p.e.PckA.kcat)*x(a.x.OAA)/(x(a.x.OAA)+p(a.p.e.PckA.Koaa)*(1+x(a.x.PEP)/p(a.p.e.PckA.Kpep)));
        
        % Pdh kinetics: MWC, substrate: PYR, inhibitors: GLX, PYR
        f(a.f.e.Pdh) = x(a.x.Pdh)*p(a.p.e.Pdh.kcat)*x(a.x.PYR)/p(a.p.e.Pdh.Kpyr)*(1+x(a.x.PYR)/p(a.p.e.Pdh.Kpyr))^(p(a.p.e.Pdh.n)-1)/((1+x(a.x.PYR)/p(a.p.e.Pdh.Kpyr))^p(a.p.e.Pdh.n)+p(a.p.e.Pdh.L)*(1+x(a.x.GLX)/p(a.p.e.Pdh.Kglx)+x(a.x.PYR)/p(a.p.e.Pdh.KpyrI))^p(a.p.e.Pdh.n));
        
        % PfkA kinetics: MWC, substrate: G6P, inhibitor: PEP
        f(a.f.e.PfkA) = x(a.x.PfkA)*p(a.p.e.PfkA.kcat)*x(a.x.G6P)/p(a.p.e.PfkA.Kg6p)*(1+x(a.x.G6P)/p(a.p.e.PfkA.Kg6p))^(p(a.p.e.PfkA.n)-1)/((1+x(a.x.G6P)/p(a.p.e.PfkA.Kg6p))^p(a.p.e.PfkA.n)+p(a.p.e.PfkA.L)*(1+x(a.x.PEP)/p(a.p.e.PfkA.Kpep))^p(a.p.e.PfkA.n));
        
        % Ppc kinetics: MWC, substrate: PEP, activator: FBP
        f(a.f.e.Ppc) = x(a.x.Ppc)*p(a.p.e.Ppc.kcat)*x(a.x.PEP)/p(a.p.e.Ppc.Kpep)*(1+x(a.x.PEP)/p(a.p.e.Ppc.Kpep))^(p(a.p.e.Ppc.n)-1)/((1+x(a.x.PEP)/p(a.p.e.Ppc.Kpep))^p(a.p.e.Ppc.n)+p(a.p.e.Ppc.L)/(1+x(a.x.FBP)/p(a.p.e.Ppc.Kfbp))^p(a.p.e.Ppc.n));
        
        % PpsA kinetics: MWC, substrate: PYR, inhibitor: PEP
        f(a.f.e.PpsA) = x(a.x.PpsA)*p(a.p.e.PpsA.kcat)*x(a.x.PYR)/p(a.p.e.PpsA.Kpyr)*(1+x(a.x.PYR)/p(a.p.e.PpsA.Kpyr))^(p(a.p.e.PpsA.n)-1)/((1+x(a.x.PYR)/p(a.p.e.PpsA.Kpyr))^p(a.p.e.PpsA.n)+p(a.p.e.PpsA.L)*(1+x(a.x.PEP)/p(a.p.e.PpsA.Kpep))^p(a.p.e.PpsA.n));
        
        % PykF kinetics: MWC, substrate: PEP, activator: FBP
        f(a.f.e.PykF) = x(a.x.PykF)*p(a.p.e.PykF.kcat)*x(a.x.PEP)/p(a.p.e.PykF.Kpep)*(1+x(a.x.PEP)/p(a.p.e.PykF.Kpep))^(p(a.p.e.PykF.n)-1)/((1+x(a.x.PEP)/p(a.p.e.PykF.Kpep))^p(a.p.e.PykF.n)+p(a.p.e.PykF.L)/(1+x(a.x.FBP)/p(a.p.e.PykF.Kfbp))^p(a.p.e.PykF.n));
               
        
        %% Gene expression rates
        
        % aceBAK expression: sum of the following three kinetics
        % MM plus basal expression, substrate: Cra
        % MM plus basal expression, substrate: Crpcamp
        % MWC-like, substrate: IclR, activator: GLX, inhibitor: PYR 
        % aceB and aceK expression are coupled to aceA expression
        % with constant factors
        f(a.f.g.aceA) = p(a.p.bm.k_expr)*mu*((1-x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.aceBAK.Kcra)))*p(a.p.g.aceBAK.vcra_unbound)...
            +x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.aceBAK.Kcra))*p(a.p.g.aceBAK.vcra_bound)...
            +(1-x(a.x.CrpcAMP)/(x(a.x.CrpcAMP)+p(a.p.g.aceBAK.Kcrp)))*p(a.p.g.aceBAK.vcrp_unbound)...
            +x(a.x.CrpcAMP)/(x(a.x.CrpcAMP)+p(a.p.g.aceBAK.Kcrp))*p(a.p.g.aceBAK.vcrp_bound)...
            +p(a.p.g.aceBAK.kcat_iclr)*x(a.x.IclR)*(1-(p(a.p.g.aceBAK.DNA)/p(a.p.g.aceBAK.KDNA))...
            *(1+x(a.x.PYR)/p(a.p.g.aceBAK.KPprime))/(1+(x(a.x.GLX)/p(a.p.g.aceBAK.KG)...
            *(1+x(a.x.GLX)/p(a.p.g.aceBAK.KG)))/p(a.p.g.aceBAK.L)...
            +p(a.p.g.aceBAK.DNA)/p(a.p.g.aceBAK.KDNA)+x(a.x.PYR)/p(a.p.g.aceBAK.KP)...
            +p(a.p.g.aceBAK.DNA)*x(a.x.PYR)/p(a.p.g.aceBAK.KDNA)/p(a.p.g.aceBAK.KPprime))));
        f(a.f.g.aceB) = p(a.p.g.aceBAK.aceBfactor)*f(a.f.g.aceA);
        f(a.f.g.aceK) = p(a.p.g.aceBAK.aceKfactor)*f(a.f.g.aceA);
        
        % acoa2act kinetics: constitutive expression
        f(a.f.g.acoa2act) = 0;
        
        % acs kinetics: Hill plus basal expression, substrate: Crpcamp
        f(a.f.g.acs) = p(a.p.bm.k_expr)*mu*((1-x(a.x.CrpcAMP)^p(a.p.g.acs.n)/(x(a.x.CrpcAMP)^p(a.p.g.acs.n)+p(a.p.g.acs.Kcrp)^p(a.p.g.acs.n)))*p(a.p.g.acs.vcrp_unbound)+x(a.x.CrpcAMP)^p(a.p.g.acs.n)/(x(a.x.CrpcAMP)^p(a.p.g.acs.n)+p(a.p.g.acs.Kcrp)^p(a.p.g.acs.n))*p(a.p.g.acs.vcrp_bound));

        % akg2mal kinetics: Hill plus basal expression, substrate: Crpcamp
        f(a.f.g.akg2mal) = p(a.p.bm.k_expr)*mu*((1-x(a.x.CrpcAMP)^p(a.p.g.akg2mal.n)/(x(a.x.CrpcAMP)^p(a.p.g.akg2mal.n)+p(a.p.g.akg2mal.Kcrp)^p(a.p.g.akg2mal.n)))*p(a.p.g.akg2mal.vcrp_unbound)+x(a.x.CrpcAMP)^p(a.p.g.akg2mal.n)/(x(a.x.CrpcAMP)^p(a.p.g.akg2mal.n)+p(a.p.g.akg2mal.Kcrp)^p(a.p.g.akg2mal.n))*p(a.p.g.akg2mal.vcrp_bound));
        
        % campdegr kinetics: constitutive expression
        f(a.f.g.campdegr) = 0;
        
        % cra kinetics: constitutive expression
        f(a.f.g.cra) = 0;
        
        % crp kinetics: constitutive expression
        f(a.f.g.crp) = 0;
        
        % cya kinetics: constitutive expression
        f(a.f.g.cya) = 0;
        
        % emp kinetics: MM plus basal expression, substrate: Cra
        f(a.f.g.emp) = p(a.p.bm.k_expr)*mu*((1-x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.emp.Kcra)))*p(a.p.g.emp.vcra_unbound)...
            +x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.emp.Kcra))*p(a.p.g.emp.vcra_bound)...
            +(1-x(a.x.CrpcAMP)/(x(a.x.CrpcAMP)+p(a.p.g.emp.Kcrp)))*p(a.p.g.emp.vcrp_unbound)...
            +x(a.x.CrpcAMP)/(x(a.x.CrpcAMP)+p(a.p.g.emp.Kcrp))*p(a.p.g.emp.vcrp_bound));
        
        % eno kinetics: MM plus basal expression, substrate: Cra
        f(a.f.g.eno) = p(a.p.bm.k_expr)*mu*((1-x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.eno.Kcra)))*p(a.p.g.eno.vcra_unbound)+x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.eno.Kcra))*p(a.p.g.eno.vcra_bound));
        
        % fdp kinetics: MM plus basal expression, substrate: Cra
        f(a.f.g.fdp) = p(a.p.bm.k_expr)*mu*((1-x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.fdp.Kcra)))*p(a.p.g.fdp.vcra_unbound)+x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.fdp.Kcra))*p(a.p.g.fdp.vcra_bound));
        
        % gltA kinetics: Hill plus basal expression, substrate: Crpcamp
        f(a.f.g.gltA) = p(a.p.bm.k_expr)*mu*((1-x(a.x.CrpcAMP)^p(a.p.g.gltA.n)/(x(a.x.CrpcAMP)^p(a.p.g.gltA.n)+p(a.p.g.gltA.Kcrp)^p(a.p.g.gltA.n)))*p(a.p.g.gltA.vcrp_unbound)+x(a.x.CrpcAMP)^p(a.p.g.gltA.n)/(x(a.x.CrpcAMP)^p(a.p.g.gltA.n)+p(a.p.g.gltA.Kcrp)^p(a.p.g.gltA.n))*p(a.p.g.gltA.vcrp_bound));
        
        % icd kinetics: MM plus basal expression, substrate: Cra
        f(a.f.g.icd) = p(a.p.bm.k_expr)*mu*((1-x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.icd.Kcra)))*p(a.p.g.icd.vcra_unbound)+x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.icd.Kcra))*p(a.p.g.icd.vcra_bound));
        
        % iclr expression: constitutive
        f(a.f.g.iclr) = 0;
        
        % mdh kinetics: MM plus basal expression, substrate: Crpcamp
        f(a.f.g.mdh) = p(a.p.bm.k_expr)*mu*((1-x(a.x.CrpcAMP)/(x(a.x.CrpcAMP)+p(a.p.g.mdh.Kcrp)))*p(a.p.g.mdh.vcrp_unbound)+x(a.x.CrpcAMP)/(x(a.x.CrpcAMP)+p(a.p.g.mdh.Kcrp))*p(a.p.g.mdh.vcrp_bound));
        
        % me kinetics: growth rate- dependent constitutive expression
        f(a.f.g.maeAB) = (mu+p(a.p.d.k_degr))*SS_MaeAB;
        
        % pckA kinetics: MM plus basal expression, substrate: Cra
        f(a.f.g.pckA) = p(a.p.bm.k_expr)*mu*((1-x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.pckA.Kcra)))*p(a.p.g.pckA.vcra_unbound)+x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.pckA.Kcra))*p(a.p.g.pckA.vcra_bound));
        
        % pdh kinetics: MM plus basal expression, substrate: PdhR
        f(a.f.g.pdh) = p(a.p.bm.k_expr)*mu*((1-x(a.x.PdhR)/(x(a.x.PdhR)+p(a.p.g.pdh.Kpdhr)))*p(a.p.g.pdh.vpdhr_unbound)+x(a.x.PdhR)/(x(a.x.PdhR)+p(a.p.g.pdh.Kpdhr))*p(a.p.g.pdh.vpdhr_bound));
        
        % pdhr expression: constitutive
        f(a.f.g.pdhr) = 0;
        
        % pfkA kinetics: MM plus basal expression, substrate: Cra
        f(a.f.g.pfkA) = p(a.p.bm.k_expr)*mu*((1-x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.pfkA.Kcra)))*p(a.p.g.pfkA.vcra_unbound)+x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.pfkA.Kcra))*p(a.p.g.pfkA.vcra_bound));
        
        % ppc kinetics: growth rate- dependent constitutive expression
        f(a.f.g.ppc) = (mu+p(a.p.d.k_degr))*SS_Ppc;
        
        % ppsA kinetics: MM plus basal expression, substrate: Cra
        f(a.f.g.ppsA) = p(a.p.bm.k_expr)*mu*((1-x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.ppsA.Kcra)))*p(a.p.g.ppsA.vcra_unbound)+x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.ppsA.Kcra))*p(a.p.g.ppsA.vcra_bound));
        
        % pykF kinetics: MM plus basal expression, substrate: Cra
        f(a.f.g.pykF) = p(a.p.bm.k_expr)*mu*((1-x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.pykF.Kcra)))*p(a.p.g.pykF.vcra_unbound)+x(a.x.Cra)/(x(a.x.Cra)+p(a.p.g.pykF.Kcra))*p(a.p.g.pykF.vcra_bound));
               
        % EIIA kinetics: constitutive expression
        f(a.f.g.EIIA) = 0;
        
        % EIICB kinetics: constitutive expression
        f(a.f.g.EIICB) = 0;
        
        
        %% Protein degradation and dilution rates
        % Constitutively produced proteins are neither produced,
        % nor degraded, nor diluted, to keep their levels constant.
        % All other proteins degrade with a constant rate,
        % and dilute with the growth rate
        
        f(a.f.d.AceA) = (mu+p(a.p.d.k_degr))*x(a.x.AceA);
        f(a.f.d.AceB) = (mu+p(a.p.d.k_degr))*x(a.x.AceB);
        f(a.f.d.AceK) = (mu+p(a.p.d.k_degr))*x(a.x.AceK);
        f(a.f.d.Acoa2act) = 0;
        f(a.f.d.Acs) = (mu+p(a.p.d.k_degr))*x(a.x.Acs);
        f(a.f.d.Akg2mal) = (mu+p(a.p.d.k_degr))*x(a.x.Akg2mal);
        f(a.f.d.CAMPdegr) = 0;
        f(a.f.d.Cra) = 0;
        f(a.f.d.CraFBP) = 0;
        f(a.f.d.Crp) = 0;
        f(a.f.d.CrpcAMP) = 0;
        f(a.f.d.Cya)  = 0;
        f(a.f.d.Emp)  = (mu+p(a.p.d.k_degr))*x(a.x.Emp);
        f(a.f.d.Eno)  = (mu+p(a.p.d.k_degr))*x(a.x.Eno);
        f(a.f.d.Fdp)  = (mu+p(a.p.d.k_degr))*x(a.x.Fdp);
        f(a.f.d.GltA) = (mu+p(a.p.d.k_degr))*x(a.x.GltA);
        f(a.f.d.Icd) = (mu+p(a.p.d.k_degr))*x(a.x.Icd);
        f(a.f.d.Icd_P) = (mu+p(a.p.d.k_degr))*x(a.x.Icd_P);
        f(a.f.d.IclR) = 0;
        f(a.f.d.Mdh)  = (mu+p(a.p.d.k_degr))*x(a.x.Mdh);
        f(a.f.d.MaeAB) = (mu+p(a.p.d.k_degr))*x(a.x.MaeAB);
        f(a.f.d.PckA) = (mu+p(a.p.d.k_degr))*x(a.x.PckA);
        f(a.f.d.Pdh)  = (mu+p(a.p.d.k_degr))*x(a.x.Pdh);
        f(a.f.d.PdhR) = 0;
        f(a.f.d.PdhRPYR) = 0;
        f(a.f.d.PfkA) = (mu+p(a.p.d.k_degr))*x(a.x.PfkA);
        f(a.f.d.Ppc)  = (mu+p(a.p.d.k_degr))*x(a.x.Ppc);
        f(a.f.d.PpsA) = (mu+p(a.p.d.k_degr))*x(a.x.PpsA);
        f(a.f.d.PykF) = (mu+p(a.p.d.k_degr))*x(a.x.PykF);
        f(a.f.d.EIIA) = 0;
        f(a.f.d.EIIA_P) = 0;
        f(a.f.d.EIICB) = 0;
        
        
        %% Metabolite dilution rates
        % Intracellular metabolites do not degrade,
        % only dilute with the growth rate.
        
        f(a.f.d.ACoA) = mu*x(a.x.ACoA);
        f(a.f.d.AKG) = mu*x(a.x.AKG);
        f(a.f.d.cAMP) = mu*x(a.x.cAMP);
        f(a.f.d.FBP) = mu*x(a.x.FBP);
        f(a.f.d.G6P) = mu*x(a.x.G6P);
        f(a.f.d.GLX) = mu*x(a.x.GLX);
        f(a.f.d.ICT) = mu*x(a.x.ICT);
        f(a.f.d.MAL) = mu*x(a.x.MAL);
        f(a.f.d.OAA) = mu*x(a.x.OAA);
        f(a.f.d.PEP) = mu*x(a.x.PEP);
        f(a.f.d.PG3) = mu*x(a.x.PG3);
        f(a.f.d.PYR) = mu*x(a.x.PYR);
       
        
        %% Environmental interaction

        % Let the cell population grow with the actual growth rate
        f(a.f.env.growth) = x(a.x.OD)*mu;
        
        % Scale glucose uptake from the environment with the
        % actual population size 
        f(a.f.env.GLCup) = p(a.p.env.uc)*p(a.p.env.M_GLC)*x(a.x.OD)*f(a.f.pts.r4);
        
        % Scale acetate uptake from the environment with the
        % actual population size
        f(a.f.env.ACTup) = p(a.p.env.uc)*p(a.p.env.M_ACT)*x(a.x.OD)*f(a.f.e.Acs);
        
        % Scale acetate excretion to the environment with the
        % actual population size for all scenarios except 1 and 4,
        % for which glucose should remain the sole carbon source -
        % in these cases, the excreted acetate is directed to nowhere.
        
        if or (scenario==1, scenario==4),
            f(a.f.env.ACTex) = 0;
        else
            f(a.f.env.ACTex) = p(a.p.env.uc)*p(a.p.env.M_ACT)*x(a.x.OD)*f(a.f.e.Acoa2act);
        end
        

        %% Return simulation result
        
        if called_by_ode15s, % If ecoliOK_core is called by ode15s ...
            vargout = S*f';  % return the state derivatives dx = S*f' ...
        else                 % otherwise ...
            vargout = f;     % return the reaction rates f
        end
         
        
    end % of nested function ecoliOK_core

end % of main function