function T = plaph_v2_tcga_summary(tcga_root)
%PLAPH_V2_TCGA_SUMMARY Collect TCGA grid, profile and subsample results.
tcga_root = char(tcga_root);
out = fullfile(tcga_root, "summary");
if ~isfolder(out)
    mkdir(out);
end
% fixed-tau grid choices
files = dir(fullfile(tcga_root, "grid", "*", "ebic_choice.csv"));
rows = cell(numel(files),1);
for f = 1:numel(files)
    rows{f} = readtable(fullfile(files(f).folder, files(f).name), "TextType","string");
end
T = vertcat(rows{:});
if ~isempty(T)
    T = sortrows(T, {'family','tau'});
    writetable(T, fullfile(out, "fixed_tau_choices.csv"));
    % gene-level supports at each fixed tau
    genes = table();
    for i = 1:height(T)
        d = char(T.output_dir(i));
        for margin = ["failure" "censoring"]
            g = dir(fullfile(d, "*_linear_" + margin + ".csv"));
            if isempty(g), continue; end
            L = readtable(fullfile(g(1).folder, g(1).name), "TextType","string");
            sel = L(is_true(L.selected) & ~is_clinical(L.feature), :);
            if isempty(sel), continue; end
            add = table(repmat(T.family(i),height(sel),1), repmat(T.tau(i),height(sel),1), ...
                repmat(margin,height(sel),1), string(sel.(1)), ...
                'VariableNames', {'family','tau','margin','feature'});
            genes = [genes; add]; %#ok<AGROW>
        end
    end
    writetable(genes, fullfile(out, "fixed_tau_supports.csv"));
end
% profile fits
pf = dir(fullfile(tcga_root, "profile", "*", "fit_row.csv"));
P = table();
for f = 1:numel(pf)
    P = [P; readtable(fullfile(pf(f).folder, pf(f).name), "TextType","string")]; %#ok<AGROW>
end
writetable(P, fullfile(out, "profile_fits.csv"));
% subsamples: selection frequencies
sf = dir(fullfile(tcga_root, "subsample", "sub_*", "fit_row.csv"));
Ssum = table();
freq = containers.Map();
for f = 1:numel(sf)
    Ssum = [Ssum; readtable(fullfile(sf(f).folder, sf(f).name), "TextType","string")]; %#ok<AGROW>
    for margin = ["failure" "censoring"]
        g = dir(fullfile(sf(f).folder, "*_linear_" + margin + ".csv"));
        if isempty(g), continue; end
        L = readtable(fullfile(g(1).folder, g(1).name), "TextType","string");
        names = string(L.feature(is_true(L.selected) & ~is_clinical(L.feature)));
        for nm = names(:)'
            key = char(margin + ":" + nm);
            if isKey(freq, key), freq(key) = freq(key) + 1; else, freq(key) = 1; end
        end
    end
end
writetable(Ssum, fullfile(out, "subsample_fits.csv"));
if freq.Count > 0
    k = string(keys(freq))'; v = cell2mat(values(freq))';
    parts = split(k, ":", 2);
    if size(parts,2) == 1, parts = reshape(parts, 1, []); end
    F = table(parts(:,1), parts(:,2), v, v/max(height(Ssum),1), ...
        'VariableNames', {'margin','feature','count','frequency'});
    F = sortrows(F, {'margin','frequency'}, {'ascend','descend'});
    writetable(F, fullfile(out, "subsample_frequencies.csv"));
end
end


function b = is_true(v)
b = ismember(lower(strtrim(string(v))), ["1" "true"]);
end


function b = is_clinical(names)
b = ismember(string(names), ["age" "male" "stage_II" "stage_III" "stage_IV"]);
end
