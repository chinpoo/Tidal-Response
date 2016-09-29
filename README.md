# Tidal-Response
% ================= Version 2.0, by CHUAN QIN =============================
% =========================================================================
% This program is used to calculate the elastic tidal response of a planet
% with 3-D elastic and density structures in any depth ranges of the manle.
% Perturbation theory is adopted for the formulation and spectral method is
% used for solving mode couplings and tidal response. The total tidal 
% response is represented by the sum of modal responses at each order of 
% perturbation: 0, 1 and 2. Each model response is determined using
% propagator matrix method, which is based on Runge-Kutta numerical scheme.
%
% ---------------------------- Implementation -----------------------------
% 1. Initialization: model setups, variables and flags
% 2. Outer loop: for each (l1,m1) in 3-D structure
%                ----- determine mode coupling hierarchy (selection rule)
%                ----- store mode coupling information (VSH expansion)
% 3. Inner loop: for each (l,m) in tidal response (from 0th to 2nd order)
%                ----- construct matrix equation (4 categories)
%                ----- apply propagator matrix method for solution
%                ----- solving tidal response for (l,m)
% 4. Output: get total response and write it to file
%
% References: 
% 1. Qin, Chuan, Shijie Zhong, and John Wahr. "Elastic tidal response of
% a laterally heterogeneous planet: a complete perturbation formulation."
% Geophysical Journal International 207.1 (2016): 89-110.
% 2. Qin, Chuan, Shijie Zhong, and John Wahr. "A perturbation method and
% its application: elastic tidal response of a laterally heterogeneous 
% planet." Geophysical Journal International 199.2 (2014): 631-647.
% =========================================================================
