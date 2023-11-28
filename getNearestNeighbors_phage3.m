function [neighbors, infprob, ph] = getNearestNeighbors_phage3(dataXY, class, phdataXY, ph_class, radius)

%parameters
% W = width of the torus
% H = heigh of the torus
%%resultant_force = zeros(size(dataXY));
len = size(dataXY,1);
periodic_boundary_type = 1;

neighbors = cell(len,1);
dist=[];
infprob = cell(len,1);
ph=cell(len,1);
dist_th = radius.^2;

for i = 1 : len

    if periodic_boundary_type == 1
        dist = [(phdataXY(:,1)-dataXY(i,1)).^2 + (phdataXY(:,2)-dataXY(i,2)).^2];
    end
    
    dist_idx =find(dist<dist_th);
    
    if isempty(dist_idx)
        neighbors{i} = [];
        infprob{i}=[];
        ph{i}=[];
    else
        infprob{i}=rand(length(dist_idx),1);
        ph{i}=ph_class(dist_idx);
    end
    
end

end