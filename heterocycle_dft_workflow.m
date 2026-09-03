clear; clc; close all;

ORCA_EXE = 'CHANGE_ME_TO_ORCA_EXECUTABLE';

projectDir = fileparts(mfilename('fullpath'));
inputDir = projectDir;
outputDir = fullfile(projectDir, 'outputs');
figureDir = fullfile(projectDir, 'figures');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

systems = ["pyrrole", "furan", "thiophene", "pyridine"];

if strcmp(ORCA_EXE, 'CHANGE_ME_TO_ORCA_EXECUTABLE')
    error('Set ORCA_EXE before running the workflow.');
end

for i = 1:numel(systems)
    name = systems(i);
    workDir = fullfile(outputDir, name);

    if ~exist(workDir, 'dir')
        mkdir(workDir);
    end

    copyfile(fullfile(inputDir, name + ".inp"), fullfile(workDir, name + ".inp"));
    copyfile(fullfile(inputDir, name + ".xyz"), fullfile(workDir, name + ".xyz"));

    previousDir = pwd;
    cleanup = onCleanup(@() cd(previousDir));
    cd(workDir);

    command = sprintf('"%s" "%s.inp" > "%s.out"', ORCA_EXE, name, name);
    status = system(command);

    if status ~= 0
        error('ORCA calculation failed for %s.', name);
    end

    clear cleanup;
end

n = numel(systems);

Molecule = strings(n,1);
HOMO_eV = nan(n,1);
LUMO_eV = nan(n,1);
Gap_eV = nan(n,1);
IP_eV = nan(n,1);
EA_eV = nan(n,1);
Hardness_eV = nan(n,1);
Electronegativity_eV = nan(n,1);
ElectronicEnergy_Eh = nan(n,1);
ImaginaryModes = nan(n,1);
MinimumFrequency_cm1 = nan(n,1);

for i = 1:n
    name = systems(i);
    result = parse_orca_output(fullfile(outputDir, name, name + ".out"));

    Molecule(i) = capitalize(name);
    HOMO_eV(i) = result.homo_eV;
    LUMO_eV(i) = result.lumo_eV;
    Gap_eV(i) = LUMO_eV(i) - HOMO_eV(i);

    IP_eV(i) = -HOMO_eV(i);
    EA_eV(i) = -LUMO_eV(i);
    Hardness_eV(i) = (IP_eV(i) - EA_eV(i)) / 2;
    Electronegativity_eV(i) = (IP_eV(i) + EA_eV(i)) / 2;

    ElectronicEnergy_Eh(i) = result.energy_Eh;
    ImaginaryModes(i) = result.imaginary_modes;
    MinimumFrequency_cm1(i) = result.minimum_frequency_cm1;
end

results = table( ...
    Molecule, HOMO_eV, LUMO_eV, Gap_eV, ...
    IP_eV, EA_eV, Hardness_eV, Electronegativity_eV, ...
    ElectronicEnergy_Eh, ImaginaryModes, MinimumFrequency_cm1);

writetable(results, fullfile(outputDir, 'electronic_descriptors.csv'));
write_summary(fullfile(outputDir, 'run_summary.txt'), results);

plot_gaps(results, fullfile(figureDir, 'homo_lumo_gaps.png'));
plot_orbitals(results, fullfile(figureDir, 'frontier_orbital_energies.png'));
plot_descriptors(results, fullfile(figureDir, 'conceptual_dft_descriptors.png'));

disp(results);


function result = parse_orca_output(filename)

text = fileread(filename);

if ~contains(text, 'ORCA TERMINATED NORMALLY')
    error('ORCA did not terminate normally: %s', filename);
end

energyTokens = regexp( ...
    text, ...
    'FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+(?:[Ee][+-]?\d+)?)', ...
    'tokens');

if isempty(energyTokens)
    error('Final energy not found in %s.', filename);
end

result.energy_Eh = str2double(energyTokens{end}{1});

sections = regexp(text, 'ORBITAL ENERGIES', 'start');

if isempty(sections)
    error('Orbital energies not found in %s.', filename);
end

block = text(sections(end):end);

rows = regexp( ...
    block, ...
    '^\s*(\d+)\s+([0-9]+\.[0-9]+)\s+(-?\d+\.\d+(?:[Ee][+-]?\d+)?)\s+(-?\d+\.\d+(?:[Ee][+-]?\d+)?)\s*$', ...
    'tokens', ...
    'lineanchors');

if isempty(rows)
    error('Could not parse orbital energies in %s.', filename);
end

occupancy = nan(numel(rows),1);
energy_eV = nan(numel(rows),1);

for j = 1:numel(rows)
    occupancy(j) = str2double(rows{j}{2});
    energy_eV(j) = str2double(rows{j}{4});
end

occupied = find(occupancy > 1e-6);
virtual = find(occupancy <= 1e-6);

if isempty(occupied) || isempty(virtual)
    error('Occupied and virtual orbitals could not be identified in %s.', filename);
end

homoRow = occupied(end);
lumoRow = virtual(find(virtual > homoRow, 1, 'first'));

if isempty(lumoRow)
    error('LUMO could not be identified in %s.', filename);
end

result.homo_eV = energy_eV(homoRow);
result.lumo_eV = energy_eV(lumoRow);

frequencySections = regexp(text, 'VIBRATIONAL FREQUENCIES', 'start');

if isempty(frequencySections)
    result.imaginary_modes = NaN;
    result.minimum_frequency_cm1 = NaN;
    return;
end

frequencyBlock = text(frequencySections(end):end);

frequencyTokens = regexp( ...
    frequencyBlock, ...
    '^\s*\d+:\s+(-?\d+\.\d+)\s+cm\*\*-1', ...
    'tokens', ...
    'lineanchors');

if isempty(frequencyTokens)
    result.imaginary_modes = NaN;
    result.minimum_frequency_cm1 = NaN;
else
    frequencies = cellfun(@(x) str2double(x{1}), frequencyTokens);
    result.imaginary_modes = sum(frequencies < -20);
    result.minimum_frequency_cm1 = min(frequencies);
end

end


function write_summary(filename, T)

fid = fopen(filename, 'w');

fprintf(fid, 'B3LYP-D3(BJ)/def2-TZVP\n');
fprintf(fid, 'TightSCF, TightOpt, DEFGRID3\n');
fprintf(fid, 'Gas phase\n\n');

for i = 1:height(T)
    fprintf(fid, '%s\n', T.Molecule(i));
    fprintf(fid, 'Electronic energy: %.12f Eh\n', T.ElectronicEnergy_Eh(i));
    fprintf(fid, 'HOMO: %.6f eV\n', T.HOMO_eV(i));
    fprintf(fid, 'LUMO: %.6f eV\n', T.LUMO_eV(i));
    fprintf(fid, 'Gap: %.6f eV\n', T.Gap_eV(i));
    fprintf(fid, 'Approximate IP: %.6f eV\n', T.IP_eV(i));
    fprintf(fid, 'Approximate EA: %.6f eV\n', T.EA_eV(i));
    fprintf(fid, 'Hardness: %.6f eV\n', T.Hardness_eV(i));
    fprintf(fid, 'Electronegativity: %.6f eV\n', T.Electronegativity_eV(i));
    fprintf(fid, 'Imaginary modes below -20 cm^-1: %g\n', T.ImaginaryModes(i));
    fprintf(fid, 'Minimum frequency: %.3f cm^-1\n\n', T.MinimumFrequency_cm1(i));
end

fclose(fid);

end


function plot_gaps(T, filename)

fig = figure('Position', [100 100 760 480]);

x = 1:height(T);
bar(x, T.Gap_eV);

xticks(x);
xticklabels(T.Molecule);
xlabel('Heterocycle');
ylabel('HOMO-LUMO gap (eV)');
title('Calculated HOMO-LUMO Gaps');
grid on;

for i = 1:height(T)
    text(i, T.Gap_eV(i) + 0.03, sprintf('%.3f', T.Gap_eV(i)), ...
        'HorizontalAlignment', 'center');
end

exportgraphics(fig, filename, 'Resolution', 300);
close(fig);

end


function plot_orbitals(T, filename)

fig = figure('Position', [100 100 860 520]);
hold on;

for i = 1:height(T)
    plot([i-0.23 i+0.23], [T.HOMO_eV(i) T.HOMO_eV(i)], 'LineWidth', 2);
    plot([i-0.23 i+0.23], [T.LUMO_eV(i) T.LUMO_eV(i)], 'LineWidth', 2);

    text(i+0.28, T.HOMO_eV(i), sprintf('HOMO %.3f', T.HOMO_eV(i)), ...
        'VerticalAlignment', 'middle');

    text(i+0.28, T.LUMO_eV(i), sprintf('LUMO %.3f', T.LUMO_eV(i)), ...
        'VerticalAlignment', 'middle');
end

xticks(1:height(T));
xticklabels(T.Molecule);
xlabel('Heterocycle');
ylabel('Orbital energy (eV)');
title('Frontier Molecular Orbital Energies');
grid on;
hold off;

exportgraphics(fig, filename, 'Resolution', 300);
close(fig);

end


function plot_descriptors(T, filename)

fig = figure('Position', [100 100 860 520]);

x = 1:height(T);
values = [T.IP_eV, T.EA_eV, T.Hardness_eV, T.Electronegativity_eV];

bar(x, values);

xticks(x);
xticklabels(T.Molecule);
xlabel('Heterocycle');
ylabel('Energy (eV)');
title('Conceptual DFT Descriptors');
legend({'IP', 'EA', '\eta', '\chi'}, 'Location', 'best');
grid on;

exportgraphics(fig, filename, 'Resolution', 300);
close(fig);

end


function value = capitalize(value)

chars = char(value);
chars(1) = upper(chars(1));
value = string(chars);

end
