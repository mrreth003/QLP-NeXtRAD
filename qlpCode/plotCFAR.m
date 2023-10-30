function plotCFAR(data, guardCells, refCells, pfa, bw, isCFARO, isRange, isFirst, prf)
persistent himage;
% Apply CFAR
cfar = applyCFAR(data, guardCells, refCells, pfa, isCFARO, isRange, isFirst);

% Calculate Axes
rAxis = (0:size(data, 2) - 1) * 299792458 / (2 * bw * 1e6);
dAxis = linspace(- (prf * 1e3) / 2, (prf * 1e3) / 2, size(data, 1));

% Plot the CFAR output
if isFirst
    figure(4);
    himage = imagesc(rAxis, dAxis, cfar);
    xlabel('Range (m)');
    ylabel('Doppler Frequency (Hz)');
    title('NeXtRAD CFAR Quick-Look');
    colormap('gray');
else
    figure(4);
    set(himage, 'CData', cfar);

    % Adjust the axis limits
    axis([0, max(rAxis), min(dAxis), max(dAxis)]);

    % Automatically set the color axis limits based on the new data
    clim auto;

    drawnow;
end
end