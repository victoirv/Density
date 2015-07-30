function plotcmds(fname,writeimgs)
%PLOTCMDS

% Prefix filenames with name of script that generated.
s = dbstack;
base = ['./figures/',s(end).name,'_',fname];

if ~exist('./figures')
	mkdir('./figures')
end

a = get(gcf,'Children');
for i = 1:length(a)
	b = get(a(i));
	if isfield(b,'Location') % The legend object.
		set(a(i),'Color','none'); % Set legend background to transparent.
	else
		set(get(a(i),'Title'),'FontSize',14)
		set(get(a(i),'XLabel'),'FontSize',14)
		set(get(a(i),'YLabel'),'FontSize',14)
	end
end

if (nargin == 2 && ~writeimgs)
	return
end

fprintf('plotcmds: Writing (using print)\n  %s.eps\n',base);
print('-depsc',sprintf('%s.eps',base))
fprintf('  Done. \n',base);

if exist('/usr/local/bin/convert')
	fprintf('plotcmds: Writing (using convert on eps)\n  %s.png\n',base);
	% MATLAB does not allow quality parameter for png.
	% Use convert if available to create a better image.
	com = sprintf('convert -quality 100 -density 150 %s.eps %s.png',base,base);
	system(com);
	fprintf('  Done.\n',base);
else
	fprintf('plotcmds: Writing (using print)\n  %s.png\n',base);
	print('-dpng',sprintf('%s.png',base))
	fprintf('  Done.\n',base);
end
