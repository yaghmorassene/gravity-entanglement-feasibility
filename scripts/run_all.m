function run_all
% Runs all demos and organizes outputs into results folders
disp('=== Running all gravity-entanglement MATLAB demos ===');

if ~exist('results/figures','dir'), mkdir('results/figures'); end
if ~exist('results/tables','dir'), mkdir('results/tables'); end
if ~exist('results/logs','dir'), mkdir('results/logs'); end

diary(fullfile('results/logs', ['run_' datestr(now,'yyyymmdd_HHMMSS') '.log']));

addpath('demos');

demos = {
    'entanglement_feasibility_Keepone', ...
    'entanglement_feasibility_Keeptwo', ...
    'entanglement_feasibility_njpKeepfour', ...
    'entanglement_feasibility_njpKeepfive', ...
    'entanglement_figs_with_errors_and_csv_Keepthree'
};

for i = 1:numel(demos)
    fprintf('\n--- Running %s ---\n', demos{i});
    try
        feval(demos{i});
    catch ME
        fprintf(2,'Error in %s: %s\n', demos{i}, ME.message);
    end
end

diary off;
disp('=== All demos complete. Results saved in results/ ===');
end
