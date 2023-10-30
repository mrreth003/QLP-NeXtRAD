function plotRDM(data, bw, isCFAR, guardCells, refCells, pfa, isCFARO, isRange, isFirst, prf)
persistent himage
rdm = fftshift(fft(data, [], 1));

if isCFAR
    plotCFAR(abs(rdm.^2), guardCells, refCells, pfa, bw, isCFARO, isRange, isFirst);
end

% Calculate Axes
rAxis = (0:size(data, 2) - 1) * 299792458 / (2 * bw * 1e6);
dAxis = linspace(- (prf * 1e3) / 2, (prf * 1e3) / 2, size(data, 1));

% Plot the RD map
if isFirst
    figure(2);
    himage = imagesc(rAxis, dAxis, 20 * log10(abs(rdm)));
    xlabel('Range (m)');
    ylabel('Doppler Frequency (Hz)');
    title('NeXtRAD R-D Map Quick-Look');
    colormap('jet');
    colorbar;
else
    figure(2);
    set(himage, 'CData', 20 * log10(abs(rdm)));

    % Adjust the axis limits
    axis([0, max(rAxis), min(dAxis), max(dAxis)]);

    % Automatically set the color axis limits based on the new data
    clim auto;

    drawnow;
end
end