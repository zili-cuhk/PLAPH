function M = plaph_v2_manifest(experiment)
%PLAPH_V2_MANIFEST Frozen task tables for the v2 simulation campaign.
%
% Each row is one scenario. Columns: scenario_index, scenario_id, n, p, q,
% target_censor_rate, tau, rho, censoring_scale, generate_copula,
% fit_copula, fit_tau (NaN = same as tau), R, seed_root, methods, xi,
% seed_scenario. Seeds: seed_root + 1000*seed_scenario + replication, where
% seed_scenario = scenario_index except where rows must share data sets:
% the working families fitted to one generating family (copula_family) and
% the kink coefficients at one sample size (xi_sensitivity). The "main" and
% derived experiments reuse the July 2026 formal seed root (25000000) so
% that old and new estimators are compared on identical data. The spline
% pilot uses a disjoint development root (26000000).

experiment = lower(string(experiment));
core = core_design();
switch experiment
    case "main"
        M = core;
        M.methods = repmat("MPL_fixed_tau|PL",height(M),1);
    case "ablation_tau0"
        % Same data as "main" (and, last row, as "highcorr"); working copula
        % fixed at independence. The highcorr row was appended on 2026-10-22,
        % so task indices 1-360 are unchanged.
        M = [core(core.tau > 0,:); highcorr_row(core)];
        M.fit_tau = zeros(height(M),1);
        M.methods = repmat("MPL_fixed_tau",height(M),1);
    case "grpsurv_export"
        % Data for grpsurv: the 13 core designs and (appended) highcorr.
        M = [core; highcorr_row(core)];
        M.methods = repmat("DataOnly",height(M),1);
    case "highcorr"
        M = highcorr_row(core);
    case "copula_family"
        base = core(core.scenario_index == 5,:);
        gen = ["Frank";"Frank";"Frank";"Frank";"Gaussian";"Gaussian"; ...
            "Clayton";"Clayton";"Gumbel";"Gumbel"];
        fit = ["Frank";"Gaussian";"Clayton";"Gumbel";"Gaussian";"Frank"; ...
            "Clayton";"Frank";"Gumbel";"Frank"];
        M = repmat(base,numel(gen),1);
        M.scenario_index = (31:40)';
        M.generate_copula = gen;
        M.fit_copula = fit;
        % Paired design: every working family is fitted to the same data
        % sets as the correctly specified fit of its generating family.
        M.seed_scenario = [31;31;31;31;35;35;37;37;39;39];
        M.methods = repmat("MPL_fixed_tau",numel(gen),1);
    case "xi_sensitivity"
        % The same 200 data sets at every xi (within-data-set paths); they
        % are independent of D5 and D8.
        xi = [0 5 10 15 20 25.6685 30 35 40 50]';
        base = core(ismember(core.scenario_index,[5 8]),:);
        M = [repmat(base(1,:),numel(xi),1); repmat(base(2,:),numel(xi),1)];
        M.xi = [xi; xi];
        M.scenario_index = (51:(50+height(M)))';
        M.seed_scenario = [repmat(51,numel(xi),1); repmat(61,numel(xi),1)];
        M.methods = repmat("MPL_fixed_tau",height(M),1);
    case "spline_pilot"
        % Development seeds only (disjoint from all reported seeds).
        M = core(ismember(core.scenario_index,[1 4 6]),:);
        M.R(:) = 20;
        M.seed_root(:) = 26000000;
        M.methods = repmat("MPL_fixed_tau",height(M),1);
    case "smoke"
        M = core(core.scenario_index == 5,:);
        M.R(:) = 2;
        M.seed_root(:) = 27000000;
        M.methods = "MPL_fixed_tau|PL";
    otherwise
        error("PLAPH:UnknownExperiment","Unknown experiment %s.",experiment);
end
unset = isnan(M.seed_scenario);
M.seed_scenario(unset) = M.scenario_index(unset);
M.scenario_id = compose("V2_%s_S%02d_n%d_p%d_cr%02d_tau%02d_%s_%s", ...
    repmat(experiment,height(M),1),M.scenario_index,M.n,M.p,round(100*M.target_censor_rate), ...
    round(100*M.tau),M.generate_copula,M.fit_copula);
end


function M = highcorr_row(core)
% D5 with autoregressive correlation 0.70 in both covariate blocks.
M = core(core.scenario_index == 5,:);
M.scenario_index = 21;
M.rho = 0.70;
M.methods = "MPL_fixed_tau|PL";
end


function T = core_design()
rows = {
    1, 500, 400, 0.30, 0.00;
    2, 500, 400, 0.30, 0.50;
    3, 500, 400, 0.30, 0.80;
    4, 500, 400, 0.50, 0.00;
    5, 500, 400, 0.50, 0.50;
    6, 500, 400, 0.50, 0.80;
    7, 700, 400, 0.50, 0.00;
    8, 700, 400, 0.50, 0.50;
    9, 700, 400, 0.50, 0.80;
   10, 500, 200, 0.50, 0.50;
   11, 500, 800, 0.50, 0.50;
   12, 500, 800, 0.50, 0.00;
   13, 500, 800, 0.50, 0.80};
T = cell2table(rows,'VariableNames', ...
    {'scenario_index','n','dimension','target_censor_rate','tau'});
T.p = T.dimension;
T.q = T.dimension;
T.rho = repmat(0.25,height(T),1);
T.censoring_scale = arrayfun(@censoring_scale,T.target_censor_rate,T.tau);
T.generate_copula = repmat("Frank",height(T),1);
T.fit_copula = repmat("Frank",height(T),1);
T.fit_tau = NaN(height(T),1);
T.R = repmat(200,height(T),1);
T.seed_root = repmat(25000000,height(T),1);
T.xi = repmat(25.6685,height(T),1);
T.seed_scenario = NaN(height(T),1);
T.methods = repmat("MPL_fixed_tau|PL",height(T),1);
T.scenario_id = strings(height(T),1);
end


function s = censoring_scale(rate,tau)
% Calibrated in July 2026 for h0T(t)=2t (failure_scale=1); unchanged.
if abs(rate-0.50) < 1e-12
    s = 1.3;
elseif abs(rate-0.30) < 1e-12 && abs(tau) < 1e-12
    s = 15.9;
else
    s = 15.7;
end
end
