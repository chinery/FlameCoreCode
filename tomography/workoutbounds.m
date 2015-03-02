% load in cam details and an simage for each camera using recon
% 
for i = size(simages,3):-1:1
    edges(:,:,i) = edge(simages(:,:,i),'log');
% edges(:,:,i) = simages(:,:,i) > 0.5;
end

% edges = zeros(size(simages));
% edges([1 size(edges,1)],[1 size(edges,2)],:) = 1;

figure(1); clf; hold on;
colours = hsv(6);
for i = 1:size(simages,3)
    P = cam(i).K * [cam(i).R cam(i).R*-cam(i).cen];
    centre = cam(i).cen;
    
    for k = 1:size(simages,2) % for each column of image
        for j = 1:size(simages,1) % for each pixel in column (row)
            pixelcount = j+(k-1)*size(simages,1)+(i-1)*size(simages,2)*size(simages,1);
            if(edges(pixelcount) > 0.1)
                xy = [(k-1)/imageRatio; (j-1)/imageRatio];
                
                pa = centre;
                apb = P\[xy; 1];
                apb(1:4) = apb(1:4)./apb(4);
                pbv(1:3) = apb(1:3)-centre;
                pbv = pbv/norm(pbv);
                if(norm(centre' + pbv) > norm(centre')) % sometimes the projection moves away from the origin. HACK THAT.
                    pbv = pbv*-1;
                end
                pb = centre' + 10*pbv;
                line = [pa';pb];
                line = (newR*[pa'; pb]')';
                
                
                
                pa = line(1,:);
                pb = line(2,:);
                plot3([pa(1,1) pb(1,1)], [pa(1,2) pb(1,2)], [pa(1,3) pb(1,3)],'Color',colours(i,:));
                plot3(pa(1,1), pa(1,2), pa(1,3),'x','Color',colours(i,:));
            end
        end
    end
    
end