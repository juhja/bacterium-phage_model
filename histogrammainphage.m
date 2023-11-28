function histogrammainphage(phage_array,CRISPR_array, no_phage, no_crispr, simsteps)

figure;
s1=subplot(2,2,[1,2]);
bar(0:no_phage, histc(sum(~isnan(phage_array),2), 0:no_phage)/size(phage_array,1));
xlabel('number of inserted phages');
ylabel('frequency');
xlim([0,no_phage]);
title(s1,['After ',num2str(simsteps), ' steps' ])
%ylim([0,0.15]);
subplot(2,2,3);
bar(0:no_crispr, histc(sum(~isnan(CRISPR_array),2), 0:no_crispr)/size(CRISPR_array,1));
xlabel('number of CRISPR spacers');
ylabel('frequency');
xlim([0,no_crispr]);
%ylim([0,0.15]);
print(['histgram',num2str(simsteps) ], '-djpeg');
end

