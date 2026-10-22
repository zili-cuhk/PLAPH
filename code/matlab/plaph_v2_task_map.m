function [row, chunk, rep_start, rep_end, total] = plaph_v2_task_map(M, task_index, chunk_size)
%PLAPH_V2_TASK_MAP Map a 1-based task index to (scenario row, chunk).
chunks_per_row = ceil(M.R/chunk_size);
cumulative = cumsum(chunks_per_row);
total = cumulative(end);
if nargin < 2 || isempty(task_index)
    row = []; chunk = []; rep_start = []; rep_end = [];
    return;
end
if task_index < 1 || task_index > total || fix(task_index) ~= task_index
    error("PLAPH:TaskIndex","Task index %d outside 1..%d.",task_index,total);
end
row = find(task_index <= cumulative, 1, "first");
previous = 0;
if row > 1
    previous = cumulative(row-1);
end
chunk = task_index - previous;
rep_start = (chunk-1)*chunk_size + 1;
rep_end = min(chunk*chunk_size, M.R(row));
end
