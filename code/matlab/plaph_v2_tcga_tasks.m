function tasks = plaph_v2_tcga_tasks(stage)
%PLAPH_V2_TCGA_TASKS Ordered task names for the TCGA LSF arrays.
%   stage 1: fixed-tau EBIC grids (must finish first)
%   stage 2: profile fits and 200 stratified subsamples
%   stage 3: additive-component and screened-dimension sensitivity at tau=0
if nargin < 1
    stage = 1;
end
if stage == 1
    two_sided = [-0.7 -0.5 -0.3 0 0.3 0.5 0.7];
    one_sided = [0.3 0.5 0.7];
    tasks = [ "grid_Frank_" + compose("%.2f",two_sided'); ...
              "grid_Gaussian_" + compose("%.2f",two_sided'); ...
              "grid_Clayton_" + compose("%.2f",one_sided'); ...
              "grid_Gumbel_" + compose("%.2f",one_sided')];
elseif stage == 3
    kappas = [0.02 0.015 0.01 0.0075 0.005 0.0025];
    tasks = ["kpath_" + compose("%.4f",kappas'); "wexp"; ...
        "dim_250_200"; "dim_500_100"];
else
    tasks = ["profile_Frank";"profile_Gaussian";"profile_Clayton";"profile_Gumbel"; ...
        "subsample_" + string((1:200)')];
end
end
