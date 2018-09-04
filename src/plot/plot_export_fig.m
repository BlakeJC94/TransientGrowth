function plot_export_fig(h,filename,pdfwidth,widthheightratio,pdflabelsize)
%plot_export_fig - Save a figure as a .png .eps and .fig and adjust font 
% sizes according to the expected graph size in the .pdf file
% Modified from Ozzy 11/08/2011
%
% Syntax: plot_export_fig(h,filename,pdfwidth,widthheightratio,pdflabelsize)
%
% Input:
% 
%    filename = file name WITHOUT EXTENSION
%    h = list of plot handles in the current figure whose marker size or
%      line width should be modified. If h is -1, then all markers/lines 
%      in the figure have their sizes changed
%    pdfwidth = graph width in the pdf file
%    widthheightratio = 16/9, 4/3, Pi whatever you prefer
%    pdflabelsize = desired label font size in the pdf file. Axis tick
%      labels and legends have a fontsize of 0.8xpdflabelsize 
% 
%  - Current figure should be the one we want to export 
%      (i.e. function acts on 'gcf')
%  - at the end of the function axis object is in manual mode. Type 'axis
%      auto' to put them in 'auto mode'. Units of the gcf are also 
%      modified to "centimeters"
%  - WARNING: Not tested for 3D plots. It might work though.
%  - Default values for A4 LaTeX document: 
%      plot_export_fig(-1, [name], 14, 7/5, 18);
% 
% 


% line width in ? (between 1 and 3)
linewidth = 0.8;  
% % size of markers in pt (1-20)
% markersize = 6; 
% ratio between label font size and axis label font size (between 0 and 1)
labeloveraxislabelratio = 0.8;  
% line width of the grid (between 1 and 3)
gridwidth = 0.5; 
% resolution in dpi (between 150 and 600)
resolution = 300;  


% Set resulotion for export
resolution_str = ['-r' num2str(resolution)];

% Work out new dimensions
dimensions = [pdfwidth (pdfwidth/widthheightratio)];  
% Work out axis label font size
labelsize = pdflabelsize;  
axeslabelsize = floor(labeloveraxislabelratio*labelsize);


hold on;
axis manual;  % to prevent scale change when resizing


% Change units of the gcf to centimeters
set(gcf, 'Units','centimeters') 
% change font to LaTeX font
set(gca,'FontName','Latin Modern Math');  
% Set line width of axis and grid (not plots)
set(gca,'linewidth',gridwidth); 
% Set display size          
set(gcf,'position',[0,0, dimensions]); 
% Prevent redimensioning by the print function
set(gcf,'PaperPositionMode','auto');  
% Set font size of the axes tick labels
set(gca,'fontsize',axeslabelsize); 
% Set axis label font size
set(get(gca,'XLabel'),'fontsize',labelsize); 
set(get(gca,'YLabel'),'fontsize',labelsize);
set(get(gca,'ZLabel'),'fontsize',labelsize); 
% Set title font size
set(get(gca,'Title'),'fontsize',labelsize);  


% if plot handles that need line width change is -1, change all
if (h == -1)                                
    h = get(gca,'Children');                
end

for i=1:length(h)
    % set line width
    set(h(i), 'linewidth' , linewidth);   
%     % set marker size   
%     set(h(i), 'MarkerSize' , markersize);    
end


hold off;


% Rescale so dont cut off xLabel
set(gca,'OuterPosition',[0.01 0.01 0.99 0.99])
% export .eps (In fact, eps does not really care about the resolution)
print('-depsc2',resolution_str, filename);  
% export .png
print('-dpng',resolution_str, filename);
% save the fig for backup     
savefig(filename); 



end