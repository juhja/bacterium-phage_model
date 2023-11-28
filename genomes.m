function genomes(CRISPR_array,phage_array,simsteps)
figure;
 ossz_mx=zeros(size(phage_array));
 for i=1:size(phage_array,1)
 ossz_mx(i,CRISPR_array(i,~isnan(CRISPR_array(i,:))))=ossz_mx(i,CRISPR_array(i,~isnan(CRISPR_array(i,:))))+1;
 ossz_mx(i,phage_array(i,~isnan(phage_array(i,:))))=ossz_mx(i,phage_array(i,~isnan(phage_array(i,:))))+2;
 end
 colormap([1 1 1; 0 1 0; 1 0 0; 0 0 1]);
 imagesc(ossz_mx, [0 3]);
 caxis([0 4]);
 hcb=colorbar('YTickLabel', {'no interaction', 'CRISPR resistance', 'temperate insertion', 'both'}, 'YTick', [0.5 1.5 2.5 3.5]);
 set(hcb, 'YTickMode', 'manual');
 ylabel('cells');
 xlabel('phages');
 title(['genomes at ' , num2str(simsteps), ' steps'])
 title(hcb, 'states');
 print(['genomes',num2str(simsteps) ], '-djpeg');


