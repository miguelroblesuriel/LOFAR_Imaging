% Replot saved map figures with x = l and y = m, without rerunning the imaging.
% Saves each new figure as .fig and .png in the folder 'replotted'.
clear all
close all
clc

%% Rebuild the l, m grid (same FracAngle as in the imaging script)
load lofar_DSP_data_1.mat
c      = 299792458;
lambda = c / freq;
dist = sqrt(  (meshgrid(poslocal(:,1)) - meshgrid(poslocal(:,1)).').^2 ...
    + (meshgrid(poslocal(:,2)) - meshgrid(poslocal(:,2)).').^2);
D  = max(dist(:));
dl = 0.25 * lambda / D;
l  = 0:dl:1;
l  = [-fliplr(l(2:end)), l]';
m  = l;
len = length(l);

%% Figures to replot (edit this list)
files = {'periodogram.fig', 'music_d2.fig', 'music_d3.fig', 'music_d20.fig', ...
    'music_d80.fig', 'music_d150.fig', 'music_d225.fig', 'music_d250.fig', ...
    'MVDR.fig', 'AAR.fig', 'AAR_modified_weights.fig'};

outDir = 'replotted';
if ~exist(outDir, 'dir'), mkdir(outDir); end

for q = 1:numel(files)
    [~, name] = fileparts(files{q});

    % --- read the data back from the saved figure -----------------------
    fig  = openfig(files{q}, 'invisible');
    hImg = findobj(fig, 'Type', 'image');
    if isempty(hImg)
        warning('%s: no image found, skipped.', files{q});
        close(fig);
        continue
    end
    hImg  = hImg(1);
    img   = hImg.CData;
    xdata = hImg.XData;
    cmap  = colormap(fig);
    clim_ = get(findobj(fig, 'Type', 'axes', '-not', 'Tag', 'Colorbar'), 'CLim');
    if iscell(clim_), clim_ = clim_{1}; end
    close(fig);

    if ~isequal(size(img), [len len])
        warning('%s: size %dx%d does not match the %dx%d grid, skipped.', ...
            files{q}, size(img,1), size(img,2), len, len);
        continue
    end

    % --- orientation ---------------------------------------------------------
    % Old figures: imagesc(img) -> XData = [1 len], rows = l, cols = m.
    % New figures: imagesc(l,m,img.') -> XData starts at l(1): already correct.
    if abs(xdata(1) - l(1)) < 1e-6
        plotImg = img;              % already stored as rows = m, cols = l
    else
        plotImg = img.';            % rows = l, cols = m  ->  transpose
    end

    % --- replot ----------------------------------------------------------------
    nf = figure('Visible', 'off');
    h  = imagesc(l, m, plotImg);
    set(h, 'AlphaData', ~isnan(plotImg));
    axis xy equal tight;
    colormap(cmap);
    caxis(clim_);
    colorbar;
    set(gca, 'XDir', 'reverse');    
    xlabel('East \leftarrow  l  \rightarrow West');
    ylabel('South \leftarrow  m  \rightarrow North');
    title(strrep(name, '_', '\_'));

    % --- save the new figure ---------------------------------------------------
    savefig(nf, fullfile(outDir, [name '_lm.fig']));
    exportgraphics(nf, fullfile(outDir, [name '_lm.png']), 'Resolution', 200);
    close(nf);
    fprintf('Saved %s_lm.fig and .png\n', name);
end