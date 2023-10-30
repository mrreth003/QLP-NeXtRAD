function plotRTI(data, bw, prf, isFirst)
% Define a handle for the image
persistent himage;

% Calculate Axes
rAxis = (0:size(data, 2) - 1) * 299792458 / (2 * bw * 1e6);
tAxis = (0:size(data, 1) - 1) * size(data, 1) / (prf * 1e3);


if isFirst
    figure(1);
    himage = imagesc(rAxis, tAxis, 20 * log10(abs(data)));
    xlabel('Range (m)');
    ylabel('Pulse Number');
    title('NeXtRAD RTI Quick-Look');
    colormap('jet');
    colorbar;
else
    figure(1);
    set(himage, 'CData', 20 * log10(abs(data)));

    % Adjust the axis limits
    axis([0, max(rAxis), 0, max(tAxis)]);

    % Automatically set the color axis limits based on the new data
    clim auto;

    drawnow;
end
end

