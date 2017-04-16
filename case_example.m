function mpc = case_example
%CASE9    Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from Joe H. Chow's book, p. 70.

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
	2	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
	3	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
	4	1	90	30	0	0	1	1	0	345	1	1.1	0.9;
	5	1	100	35	0	0	1	1	0	345	1	1.1	0.9;
	6	1	125	50	0	0	1	1	0	345	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	300	-300	1	100	1	80	10	0	0	0	0	0	0	0	0	0	0	0;
	2	150	0	300	-300	1	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
	3	75	0	300	-300	1	100	1	100	10	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	4	0.015    0.01	0   70	70	70	0	0	1	-360	360;
	4	5	0.015	 0.01	0   0	0	0	0	0	1	-360	360;
	3	5	0.015    0.01	0   0	0	0	0	0	1	-360	360;
	5	6	0.015	 0.01	0   0	0	0	0	0	1	-360	360;
	6	2	0.015    0.01	0   0	0	0	0	0	1	-360	360;
	6	4	0.015	 0.01	0   0	0	0	0	0	1	-360	360;
];

% rate is flow limit

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	1500	0	3	0.11	5	150;
	2	2000	0	3	0.085	1.2	600;
	2	3000	0	3	0.1225	1	335;
];
