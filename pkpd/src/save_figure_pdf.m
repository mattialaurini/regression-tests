function save_figure_pdf(file)
%SAVE_FIGURE_PDF Save the current figure (vector PDF when the extension is .pdf).
if exist('exportgraphics', 'file') == 2 || exist('exportgraphics', 'builtin') == 5
    exportgraphics(gcf, file, 'ContentType', 'vector');
else
    [~, ~, ext] = fileparts(file);
    if strcmpi(ext, '.pdf'), print(gcf, file, '-dpdf'); else, saveas(gcf, file); end
end
end
