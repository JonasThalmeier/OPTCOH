function scatterplot_func(X_distorted_AWGN,Y_distorted_AWGN,X_CD_rec,Y_CD_rec,X_freq_rec,Y_freq_rec,X_eq,Y_eq,X_sync,Y_sync,scatter_vec,prjcname,savefigure,M);
if ~exist(prjcname, 'dir')
    mkdir(prjcname);
end
num_bins = 50; % Number of bins for histogram
num_points = 1e4; % Number of points for scatter plot
jumps = round(0.5 * length(X_eq) / num_points); % Step size for plotting
axlim = sqrt(M) + 1; % Axis limits for scatter plot
pointsize = 10; % Size of points in scatter plot
if scatter_vec(1) == 1
    figure;
    subplot(1,2,1)
    [density, ~, ~, binX, binY] = histcounts2(real(X_distorted_AWGN(1:16*jumps:end)), imag(X_distorted_AWGN(1:16*jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_distorted_AWGN(1:16*jumps:end)), imag(X_distorted_AWGN(1:16*jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('X-Polarization');
    xlabel('Re');
    ylabel('Im');
    grid on;
    axis square;
    hold off;
    subplot(1,2,2)
    [density, ~, ~, binX, binY] = histcounts2(real(Y_distorted_AWGN(1:16*jumps:end)), imag(Y_distorted_AWGN(1:16*jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(Y_distorted_AWGN(1:16*jumps:end)), imag(Y_distorted_AWGN(1:16*jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('Y-Polarization');
    xlabel('Re');
    ylabel('Im');
    grid on;
    axis square;
    hold off;
    sgtitle('Receiver input');

    if savefigure == 1
        savefig(fullfile(prjcname,'OPTCOH_scatter_plot_RX_in'));
    end
end
if scatter_vec(2) == 1
    figure;
    subplot(1,2,1)
    [density, ~, ~, binX, binY] = histcounts2(real(X_CD_rec(1:4*jumps:end)), imag(X_CD_rec(1:4*jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_CD_rec(1:4*jumps:end)), imag(X_CD_rec(1:4*jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('X-Polarization');
    xlabel('Re');
    ylabel('Im');
    grid on;
    axis square;
    hold off;
    subplot(1,2,2)
    [density, ~, ~, binX, binY] = histcounts2(real(Y_CD_rec(1:4*jumps:end)), imag(Y_CD_rec(1:4*jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(Y_CD_rec(1:4*jumps:end)), imag(Y_CD_rec(1:4*jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('Y-Polarization');
    xlabel('Re');
    ylabel('Im');
    grid on;
    axis square;
    hold off;
    sgtitle('Chromatic Dispersion compensation output');
    if savefigure == 1
        savefig(fullfile(prjcname,'OPTCOH_scatter_plot_CD_comp_out'));
    end
end
if scatter_vec(3) == 1
    figure;
    subplot(1,2,1)
    [density, ~, ~, binX, binY] = histcounts2(real(X_freq_rec(round(end/2):jumps:end)), imag(X_freq_rec(round(end/2):jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_freq_rec(round(end/2):jumps:end)), imag(X_freq_rec(round(end/2):jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('X-Polarization');
    xlabel('Re');
    ylabel('Im');
    grid on;
    axis square;
    hold off;
    subplot(1,2,2)
    [density, ~, ~, binX, binY] = histcounts2(real(Y_freq_rec(round(end/2):jumps:end)), imag(Y_freq_rec(round(end/2):jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(Y_freq_rec(round(end/2):jumps:end)), imag(Y_freq_rec(round(end/2):jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('Y-Polarization');
    xlabel('Re');
    ylabel('Im');
    grid on;
    axis square;
    hold off;
    sgtitle('Frequency compensation output');
    if savefigure == 1
        savefig(fullfile(prjcname,'OPTCOH_scatter_plot_freq_comp_out'));
    end
end
if scatter_vec(4) == 1
    figure;
    subplot(1,2,1)
    [density, ~, ~, binX, binY] = histcounts2(real(X_eq(round(end/2):jumps:end)), imag(X_eq(round(end/2):jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_eq(round(end/2):jumps:end)), imag(X_eq(round(end/2):jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('X-Polarization');
    xlabel('Re');
    ylabel('Im');
    grid on;
    axis square;
    hold off;
    subplot(1,2,2)
    [density, ~, ~, binX, binY] = histcounts2(real(Y_eq(round(end/2):jumps:end)), imag(Y_eq(round(end/2):jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(Y_eq(round(end/2):jumps:end)), imag(Y_eq(round(end/2):jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('Y-Polarization');
    xlabel('Re');
    ylabel('Im');
    grid on;
    axis square;
    hold off;
    sgtitle('Equalizer output');
    if savefigure == 1
        savefig(fullfile(prjcname,'OPTCOH_scatter_plot_EQ_out'));
    end
end
if scatter_vec(5) == 1
    figure;
    subplot(1,2,1)
    [density, ~, ~, binX, binY] = histcounts2(real(X_sync(round(end/2):jumps:end)), imag(X_sync(round(end/2):jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_sync(round(end/2):jumps:end)), imag(X_sync(round(end/2):jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('X-Polarization');
    xlabel('Re');
    ylabel('Im');
    grid on;
    axis square;
    hold off;
    subplot(1,2,2)
    [density, ~, ~, binX, binY] = histcounts2(real(Y_sync(round(end/2):jumps:end)), imag(Y_sync(round(end/2):jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(Y_sync(round(end/2):jumps:end)), imag(Y_sync(round(end/2):jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('Y-Polarization');
    xlabel('Re');
    ylabel('Im');
    grid on;
    axis square;
    hold off;
    sgtitle('Phase compensation output');

    if savefigure == 1
        savefig(fullfile(prjcname,'OPTCOH_scatter_plot_phase_comp_out'));
    end
end
end

