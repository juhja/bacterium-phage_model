function Matrix(CRISPR_array,bact_class,phage_array,SPECIES_NUM,no_phage,simsteps)
        
    figure;
    
    subplot(1,2,1)
    temp_mx=zeros(SPECIES_NUM,no_phage);
    for i=1:SPECIES_NUM
       for j=1:no_phage
           temp_mx(i,j)=sum(sum(phage_array(bact_class==i,:)==j))/sum(bact_class==i); 
       end
    end
    imagesc(temp_mx);
    title([num2str(simsteps), ' steps; temperate phage ratios']);%numbers');
    xlabel('phage type');
    ylabel('cell type');
    colorbar;
    
    subplot(1,2,2)
    CRISPR_mx=zeros(SPECIES_NUM,no_phage);
    for i=1:SPECIES_NUM
       for j=1:no_phage
           CRISPR_mx(i,j)=sum(sum(CRISPR_array(bact_class==i,:)==j))/sum(bact_class==i); 
       end
    end
    imagesc(CRISPR_mx);
    title('CRISPR spacer ratios');%numbers');
    xlabel('phage type');
    ylabel('cell type');
    colorbar;
    print(['matrix',num2str(simsteps) ], '-djpeg');

end
