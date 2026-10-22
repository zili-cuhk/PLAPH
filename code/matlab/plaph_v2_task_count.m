function total = plaph_v2_task_count(experiment, chunk_size)
%PLAPH_V2_TASK_COUNT Number of LSF array tasks for an experiment.
if nargin < 2
    chunk_size = 5;
end
M = plaph_v2_manifest(experiment);
[~,~,~,~,total] = plaph_v2_task_map(M, [], chunk_size);
fprintf("%s: %d tasks (chunk size %d)\n", string(experiment), total, chunk_size);
end
