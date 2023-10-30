function result = applyCFAR(data, guardCellsRange, refCellsRange, pfa, isCFARO, isRange, isFirst) % Sliding window in columns
    [numRows, numCols] = size(data);
    result = zeros(size(data));

    if isRange
        cfar = data;
        thresh = zeros(size(data));
        for i = 1:numRows
            for j = 1:numCols
                % Adjust reference cell extraction based on the position in the matrix
                if j < guardCellsRange + refCellsRange + 1
                    % Handle left edge cases with adjusted reference cell extraction
                    refCells = cfar(i, j + 1 + guardCellsRange : j + refCellsRange + guardCellsRange);
                elseif j > numCols - guardCellsRange - refCellsRange
                    % Handle right edge cases with adjusted reference cell extraction
                    refCells = cfar(i, j - refCellsRange - guardCellsRange : j - 1 - guardCellsRange);
                else
                    % Standard reference cell extraction
                    refCells = [cfar(i, j - refCellsRange - guardCellsRange : j - guardCellsRange - 1), cfar(i, j + guardCellsRange + 1: j + refCellsRange + guardCellsRange)];
                end

                % CA-CFAR
                alphaCA = size(refCells, 2) * (pfa ^ (-1 / size(refCells, 2)) - 1);
                gCA = (1 / size(refCells, 2)) * sum(refCells);
                threshold = alphaCA * gCA;

                if cfar(i, j) > threshold
                    result(i, j) = 1;
                end

                thresh(i, j) = threshold;
            end
        end
    else
        cfar = data;
        thresh = zeros(size(data));
        for i = 1:numCols
            for j = 1:numRows
                % Adjust reference cell extraction based on the position in the matrix
                if j < guardCellsRange + refCellsRange + 1
                    % Handle left edge cases with adjusted reference cell extraction
                    refCells = cfar(j + 1 + guardCellsRange : j + refCellsRange + guardCellsRange, i);
                elseif j > numRows - guardCellsRange - refCellsRange
                    % Handle right edge cases with adjusted reference cell extraction
                    refCells = cfar(j - refCellsRange - guardCellsRange : j - 1 - guardCellsRange, i);
                else
                    % Standard reference cell extraction
                    refCells = [cfar(j - refCellsRange - guardCellsRange : j - guardCellsRange - 1, i); cfar(j + guardCellsRange + 1: j + refCellsRange + guardCellsRange, i)];
                end
                    
                % CA-CFAR
                alphaCA = size(refCells, 1) * (pfa ^ (-1 / size(refCells, 1)) - 1);
                gCA = (1 / size(refCells, 1)) * sum(refCells);
                threshold = alphaCA * gCA;

                if cfar(j, i) > threshold
                    result(j, i) = 1;
                end

                thresh(j, i) = threshold;
            end
        end
    end

    if isCFARO
        thresh = reshape(thresh, 1, []);
        cfar = reshape(cfar, 1, []);

        if isFirst
            figure(5);
            plot(20*log10(thresh));
            xlabel('Range (m)');
            ylabel('Power (dB)');
            title('NeXtRAD CFAR Quick-Look');
            hold on
            plot(20*log10(cfar));
            hold off
        else
            figure(5);
            clf;
            plot(20*log10(thresh));
            xlabel('Range (m)');
            ylabel('Power (dB)');
            title('NeXtRAD CFAR Quick-Look');
            hold on
            plot(20*log10(cfar));
            hold off
        end

    end
end