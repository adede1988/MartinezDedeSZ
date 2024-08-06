function [p_bin] = addRedOutline(p, alpha, colVal)
hold on 
p_bin = zeros(size(p)); 
p_bin(p<alpha) = 1; %1 = cluster!

out = p_bin; 
phor = diff(p_bin);
pver = diff(p_bin'); 
phor(phor ~= 0) = 1; 
pver(pver ~= 0) = 1; 
phor = [zeros(size(phor,2),1)'; phor]; 
pver = [zeros(size(pver,2),1)'; pver]; 
p_bin = phor + pver'; 
p_bin(p_bin>1) = 1; 


for ii = 1:size(p_bin,2)
    cur = p_bin(:,ii); 
    pnt = find(cur==1);
    if length(pnt)>0

        scatter(pnt-.5, ii-.5, 65, 'filled', 'square', ...
            'MarkerEdgeColor', colVal,...
            'MarkerFaceColor', colVal)

%     breaks = find(diff(pnt)>1); 
%     breaks = [0; breaks; length(pnt)]; 
%     for bb = 1:length(breaks)-1
%         x = [breaks(bb)+1:breaks(bb+1)]; 
%         y = ones(length(x),1)*(ii-.5);
%         x = [x, flip(x), x(1)]; 
%         y = [y; y+1; y(1)]; 
%         plot(pnt(x)-.5, y, 'color', colVal, 'linewidth', 2)
% 
% %     scatter(pnt, ii, 5, colVal, 'filled', 'square')
%     end
    end
end