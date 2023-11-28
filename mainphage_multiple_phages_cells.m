% INITIALIZAITON AND PARAMETERS
tstart=tic;

inputname='mainphage parameters.xlsx'; % input data
inputsheet='Scenario2_b5_p5'; %Scenario1,Scenario2,Scenario3,Scenario4, _b1_p60, _b5_p5,
inputsheet_phagespec='nested';%'modular';%'one-to-one';%0;% 
community_development_step=1000; % supplementary plot intervals

% read in parameters:
PLOT=xlsread(inputname, inputsheet, 'B2:B2'); % 1=plot 0=do not plot 
bacnb=xlsread(inputname,inputsheet,'B3:B3'); % number of bacteria at the beginning
SPECIES_NUM=xlsread(inputname,inputsheet,'B4:B4'); % number of cell types at the beginning
growth_step=xlsread(inputname,inputsheet,'B5:B5'); %cell division in every X steps
growth_perc=xlsread(inputname,inputsheet,'B6:B6'); % dividing ratio of the cells
death_step =xlsread(inputname,inputsheet,'B7:B7'); % cell death in every X steps
death_rate=xlsread(inputname,inputsheet,'B8:B8'); % dieing ratio of the cells
maxsimsteps=xlsread(inputname,inputsheet,'B9:B9'); % maximal iteration number
STEP_LENGTH=xlsread(inputname,inputsheet,'B10:B10'); % length of bacterial step
ground_size=xlsread(inputname,inputsheet,'B11:B11'); % this is the size of the plate in the hard_boundary=1 case
center1=xlsread(inputname,inputsheet,'B12:C12'); % middle of the plate, initial agent positions 
hard_boundary=xlsread(inputname,inputsheet,'B13:B13'); % 1 means the bacteria is placed on a plate with hard boundary. In this case the ground is finite.
no_phage=xlsread(inputname,inputsheet,'B14:B14'); % number of phage types
no_crispr=xlsread(inputname,inputsheet,'B15:B15'); % number of CRISPR spacers
scale_ph=xlsread(inputname,inputsheet,'B16:B16'); % phage colos scale for the plotting
scale_crispr=xlsread(inputname,inputsheet,'B17:B17'); % CRISPR circle size for the plotting
restmod_rate=xlsread(inputname,inputsheet,'B18:B18'); % probability for phage degradation by the restriction modification system
crisprinterg_rate=xlsread(inputname,inputsheet,'B19:B19'); % probability of integration into CRISPR
virulent_rate=xlsread(inputname,inputsheet,'B20:B20'); % probability of the phage is virulent
phageinit_kezd=xlsread(inputname,inputsheet,'B21:B21'); % number of phages at the beginning
phagenbinit=xlsread(inputname,inputsheet,'B22:B22'); % phage replication rate: X phage from 1 in the lytic cycle
phradius=xlsread(inputname,inputsheet,'B23:B23'); % radius of phage spread after replication
R=xlsread(inputname,inputsheet,'B24:B24'); % in this radius the cell can interact with the phage
phinf_rate=xlsread(inputname,inputsheet,'B25:B25'); % probability of phage infection if it is close enough to the bacteria: phinf_rate
ph_decay=xlsread(inputname,inputsheet,'B26:B26'); % a phage can work till X steps than degrades (disappears)

if inputsheet_phagespec~=0 %different phage_specificity profiles: rows are the bacteria types, columns are the phage types (if it is not uniform defined by phinf_rate)
    phinf_rate=xlsread(inputname,inputsheet_phagespec);
end

bacteria_population_envf =sprintf('phagesystem_%d_%d_%G.mat', no_phage, no_crispr, R);%name of the output

close all;

% simluation parameters:
bacnb = bacnb*SPECIES_NUM;  %skip
bact_temp_state = []; 
growth_cnt =0; %%cell division counter
death_cnt = 0; %%cell death counter
changes = 0;
changes_cnt = 0;

% baceria:
POPULATION_SIZE = ground_size*200; % maximal cell population size
bact_prop=[repmat(STEP_LENGTH', bacnb/length(STEP_LENGTH), 1)];
radius = ground_size;%/2; %%put cells in an area with that radius

% PLACE BACTERIA ON THE GROUND:
    alpha = rand(bacnb,1).*2*pi;
    random_vector = [sin(alpha), cos(alpha)];
    hossz= rand(bacnb,1)+rand(bacnb,1);
    hosszu=find(hossz>1);
    hossz(hosszu)=2-hossz(hosszu);
    positions =random_vector.*[hossz hossz].*radius+repmat(center1,bacnb,1);% rows: cells; columns: X, Y coordinates

bact_temp_state = ones(bacnb,1);
bact_class=randi([1, SPECIES_NUM],[bacnb,1]); % set bacterial types, 

bact_id=1:bacnb; % bacteria IDs

% phages:
    phage_array=nan(bacnb, no_phage);
         szorzo=600/scale_crispr^2; % multiplication factor for the plotting
    CRISPR_array=nan(bacnb, no_crispr);
    free_freq=zeros(maxsimsteps, no_phage);
    temp_freq=zeros(maxsimsteps, no_phage);
    crispr_freq=zeros(maxsimsteps, no_phage);


% phage positions:
    alpha = rand(phageinit_kezd,1).*2*pi;
    random_vector = [sin(alpha), cos(alpha)];
    hossz= rand(phageinit_kezd,1)+rand(phageinit_kezd,1);
    hosszu=find(hossz>1);
    hossz(hosszu)=2-hossz(hosszu);
    ph_positions =random_vector.*[hossz hossz].*radius+repmat(center1,phageinit_kezd,1);% rows: phages; columns: X, Y coordinates
    ph_class=randi([1, no_phage],[phageinit_kezd,1]);
    ph_age=ones(phageinit_kezd,1);


% output matrices initialisation:
popnum=nan(maxsimsteps,8+SPECIES_NUM+no_phage);
uptake_crispr=0;
uptake_temp=0;
killed_track = zeros(maxsimsteps,1);
dead_bacteria=0;

%%
% INIT THE FIGURE WINDOW:
if PLOT == 1
    hold on;
    findfigs;
    mainaxes = gca;
    set(mainaxes, 'ylim', [-8.5,4.5]);
    set(mainaxes, 'xlim', [-4.5,4.5]);
   colormap([linspace(0,1,no_phage+1)', linspace(1,0,no_phage+1)', 0*ones(no_phage+1,1)]);
   color_type = distinguishable_colors(no_phage);

       scatter(mainaxes,linspace(-4,4,scale_crispr), -8*ones(1,scale_crispr), (1:scale_crispr).^2*szorzo, linspace(scale_ph,0,scale_crispr), 'filled'); %plot CRISPR circles
       for i=1:no_phage % initial plot the phages
       plot(mainaxes,ph_positions(ph_class==i,1),ph_positions(ph_class==i,2), 'Color', color_type(i,:), 'MarkerFaceColor', color_type(i,:), 'marker', 'o','MarkerSize', 0.001, 'linestyle','none');
       end
       
       for i=1:SPECIES_NUM % change if diff number of phage and cell types
       plot(mainaxes,ph_positions(ph_class==i,1),ph_positions(ph_class==i,2),'.', 'Color', color_type(i,:), 'MarkerFaceColor', color_type(i,:), 'marker', 'o','MarkerSize', 3,'linestyle','none'); % plot initial bacteria    
       end
       
       xlabel(['circle size is proportional to the number of CRISPR spacers (0-', num2str(scale_crispr-1), ')']);
       set(gca, 'XTick', [], 'YTick', []);
       hcb1=colorbar;
       title(hcb1, 'number of integrated temperate phages');
       
    % SUBPLOT FOR free phage frequencies:
    subaxes_phage= axes('Position', [.2, .19, .25, .07]);
    set(subaxes_phage, 'ylim', [0,1]);
    set(subaxes_phage, 'xlim', [1,no_phage]);
    set(subaxes_phage, 'Fontsize', 7);
    set(get(subaxes_phage,'XLabel'),'String','phages', 'Fontsize', 7);
    set(get(subaxes_phage,'YLabel'),'String','ferquency', 'Fontsize', 7);
    
    % SUBPLOT FOR temperate phage frequencies:
    subaxes_phagetemp= axes('Position', [.2, .31, .25, .07]);
    set(subaxes_phagetemp, 'ylim', [0,1]);
    set(subaxes_phagetemp, 'xlim', [1,no_phage]);
    set(subaxes_phagetemp, 'Fontsize', 7);
    set(get(subaxes_phagetemp,'XLabel'),'String','temperate phages', 'Fontsize', 7);
    set(get(subaxes_phagetemp,'YLabel'),'String','ferquency', 'Fontsize', 7);
    
    % SUBPLOT FOR CRISP spacer frequencies:
    subaxes_phagecrispr= axes('Position', [.5, .31, .25, .07]);
    set(subaxes_phagecrispr, 'ylim', [0,1]);
    set(subaxes_phagecrispr, 'xlim', [1,no_phage]);
    set(subaxes_phagecrispr, 'Fontsize', 7);
    set(get(subaxes_phagecrispr,'XLabel'),'String','CRISPR spacers', 'Fontsize', 7);
    set(get(subaxes_phagecrispr,'YLabel'),'String','ferquency', 'Fontsize', 7);
         
end

%%
% iteration of the simulation:
for simsteps = 1:maxsimsteps
    simsteps;
    % Rule 1: each bacteria makes a random step
    % random step:
    alpha = rand(bacnb,1) * 2*pi;
    random_vector = [sin(alpha), cos(alpha)];
    positions=positions+random_vector.*[bact_prop(:,1), bact_prop(:,1)];
   
    % If the bacteria are in a Petri-dish, don't let the bacteria walk out the Petri-dish
    if hard_boundary == 1 % stays at Petri-dish
        pos = positions - repmat([0,0], bacnb,1); 
        idx = find(sum(pos.^2,2)>(ground_size.^2));
        if ~isempty(idx)
            %{
      x
x = ----- * r
     ||x||
normalisation: project the x vector to circle with radius 'r'.
            %}
            positions(idx,:) = repmat([0,0],length(idx),1)+pos(idx,:)./repmat(sqrt(sum(pos(idx,:).^2,2)),1,2)*ground_size;
        end    
    end 
    
    % Rule 2: Modelling the infection for each bacteria
    survive = ones(bacnb,1); % store survivor cells here (1 survivor, 0 dead)
     % identify the bacteria whitin an R vicinity (neighbours)
     % R=0.1; %cell can uptake the phage from that distance
     % phinfrate=0.99; %phage infection rate
    [phneighbors, infprob, ph] = getNearestNeighbors_phage3(positions, bact_temp_state, ph_positions, ph_class, R); % phages that are close enough to a bacterium 
    
    for i = 1 : bacnb
        phneighs = infprob{i}; % infection probability
        ph2=ph{i}; % phage types
        row_sum=[];
        
        for i_ph=1:no_phage
            if numel(phinf_rate)==1
                [row_act, col_act]=find(phneighs<phinf_rate & ph2==i_ph); %who can infect
            else
                [row_act, col_act]=find(phneighs<phinf_rate(bact_class(i), i_ph) & ph2==i_ph); %who can infect
            end
        row_sum=[row_sum; row_act];
        end
        
           if ~isempty(row_sum) % if there is any phage that infect
           ph_class_akt=ph2(row_sum);

           for ff=1:length(ph_class_akt) % phages that try to infect
              if sum(phage_array(i,:)==ph_class_akt(ff))==0 && sum(CRISPR_array(i,:)==ph_class_akt(ff))==0 % if bact is not resistant to the chosen phage
                
                  if rand<restmod_rate % if the cell can degrade the phage
                    if rand<crisprinterg_rate % build in CRISPR
                        CRISPR_array(i,:)=[ph_class_akt(ff),CRISPR_array(i,1:no_crispr-1)]; % new CRISPR spacer
                        uptake_crispr=uptake_crispr+1;
                    end
                    
                  else % cell cannot degrade the phage
                    if rand<virulent_rate % virulent phage
                        survive(i) = 0; % cell dies
                        % new phages:
                        ph_positions =[ph_positions; randn(phagenbinit,2)*phradius+repmat(positions(i,:),phagenbinit,1)];
                        ph_class=[ph_class; ones(phagenbinit,1)*ph_class_akt(ff)];
                        ph_age=[ph_age; zeros(phagenbinit,1)];
                    else % temperate phage
                        bact_temp_state(i)=2; 
                        phage_array(i,:)=[ph_class_akt(ff),phage_array(i,1:no_phage-1)]; % phage intergration
                        uptake_temp=uptake_temp+1;
                    end
                    
                  end
                  
              end
           end
           
           end
  
    end
       
    % delete dead bacteria:
    idx = find(survive == 1); % find the survivors, store them
    % update the variables of bacteria:
    positions = positions(idx,:); 
    bact_temp_state = bact_temp_state(idx);
    bact_class=bact_class(idx);
    phage_array=phage_array(idx,:);
    CRISPR_array=CRISPR_array(idx,:);
    bact_prop=bact_prop(idx,:);
    changes = changes + (bacnb-length(idx));
    bact_id=bact_id(idx);
    bacnb = size(positions,1);% number of bacteria
    killed_track(simsteps) = sum(survive==0);
    phagenb = size(ph_positions,1); % number of phages
    
    if bacnb==0 % end the simulation is there are no cells left
        break;
    end
    
    % store the variables of the phages:
    ph_freq_free=histc(ph_class, (1:no_phage))/length(ph_class);
    ph_freq_temp=sum(histc(phage_array, (1:no_phage),2))/bacnb;
    ph_freq_crispr=sum(histc(CRISPR_array, (1:no_phage),2))/bacnb;
    free_freq(simsteps,:)=ph_freq_free;
    temp_freq(simsteps,:)=ph_freq_temp;
    crispr_freq(simsteps,:)=ph_freq_crispr;
    
    growth_cnt = growth_cnt + 1; % division counter
    death_cnt=death_cnt+1; % death counter
    changes_cnt = changes_cnt + 1; % there was phage-bact interaction in this step
    
   % plotting section (like at the beginning): 
   if PLOT == 1 % if ther is plotting

    cla(mainaxes);
    scatter(mainaxes,linspace(-4,4,scale_crispr), -8*ones(1,scale_crispr), (1:scale_crispr).^2*szorzo, linspace(scale_ph,0,scale_crispr), 'filled');
    
      for i=1:no_phage   
       plot(mainaxes,ph_positions(ph_class==i,1),ph_positions(ph_class==i,2), 'Color', color_type(i,:), 'MarkerFaceColor', color_type(i,:), 'marker', 'o','MarkerSize', 0.001, 'linestyle','none');
      end
      
      for i=1:SPECIES_NUM % change if diff number of phage and cell types  
      plot(mainaxes,positions(bact_class==i,1),positions(bact_class==i,2),'.', 'Color', color_type(i,:), 'MarkerFaceColor', color_type(i,:), 'marker', 'o','MarkerSize', 3,'linestyle','none'); % plot bacteria
      end      
       
      bar(subaxes_phage, (1:no_phage), ph_freq_free, 'b', 'EdgeColor', 'b');
      set(subaxes_phage, 'ylim', [0,1], 'xlim', [0,no_phage], 'Fontsize', 7);
      set(get(subaxes_phage,'XLabel'),'String','phages', 'Fontsize', 7);
      set(get(subaxes_phage,'YLabel'),'String','ferquency', 'Fontsize', 7);
      
      bar(subaxes_phagetemp, (1:no_phage), ph_freq_temp, 'r', 'EdgeColor', 'r');
      set(subaxes_phagetemp, 'ylim', [0,1], 'xlim', [0,no_phage], 'Fontsize', 7);
      set(get(subaxes_phagetemp,'XLabel'),'String','temperate phages', 'Fontsize', 7);
      set(get(subaxes_phagetemp,'YLabel'),'String','ferquency', 'Fontsize', 7);
      
      bar(subaxes_phagecrispr, (1:no_phage), ph_freq_crispr, 'g', 'EdgeColor', 'g');
      set(subaxes_phagecrispr, 'ylim', [0,1], 'xlim', [0,no_phage], 'Fontsize', 7);
      set(get(subaxes_phagecrispr,'XLabel'),'String','CRISPR spacers', 'Fontsize', 7);
      set(get(subaxes_phagecrispr,'YLabel'),'String','ferquency', 'Fontsize', 7);
      
      drawnow; % update the plot
 
   end
    
%   Rule 3: division of the cells
    if growth_cnt == growth_step % if there is division step
        if bacnb < POPULATION_SIZE
            
            idx = randperm(bacnb); % randomise cell order
            growth_num = round(length(idx)*growth_perc*(1-bacnb/POPULATION_SIZE)); % sigmoid growth
                                  if growth_num == 0 && bacnb < 10  % minimal growth in csase of small cell numbers
                                  growth_num = 1;
                                  end
              
            idx = idx(1:growth_num);
            % copy the data of the new cells at the end of the cell matrices
            bacnb = bacnb+length(idx);
            bact_temp_state  = [bact_temp_state;bact_temp_state(idx)];
            bact_class=[bact_class;bact_class(idx)];
            positions = [positions;(positions(idx,:))];
            bact_id=[bact_id, bact_id(idx)];
            bact_prop=[bact_prop; bact_prop(idx, :)];
            phage_array=[phage_array; phage_array(idx,:)];
            CRISPR_array=[CRISPR_array; CRISPR_array(idx,:)];
            
        end
        growth_cnt = 0;
    end
    
            % Rule 4: natural death 
                if death_cnt == death_step  % natural death in death steps
                    if bacnb > 0
                        % randomise cell order, choose the cell which will die
                        judgement = rand(bacnb,1); 
                        life = judgement;
                        [dummy,idx2] = sort(life, 'ascend');
                        idx = idx2(1:round(bacnb*(1-death_rate)));    % sigmoid survivor 
                        dead=idx2(round(bacnb*(1-death_rate))+1:end);
                        
                        % temperate phage version 1: activation for natural cell death  
                        phnew_id=dead(bact_temp_state(dead)==2); % phages of dying cells
                         phnew_num=length(phnew_id);
                         if phnew_num~=0
                             % positions and attributes of the new phages:
                             posx=reshape(repmat(positions(phnew_id,1),[1,phagenbinit])', [phagenbinit*phnew_num, 1]);
                             posy=reshape(repmat(positions(phnew_id,2),[1,phagenbinit])', [phagenbinit*phnew_num, 1]);
                             ph_positions =[ph_positions; randn(phnew_num*phagenbinit,2)*phradius+[posx,posy]];
                             ph_age=[ph_age; zeros(phnew_num*phagenbinit,1)];
                             % choose random a temperate phage and activate it:   
                             phage_aktiv=phage_array(phnew_id,:);
                             ph_num1=sum(~isnan(phage_aktiv),2);
                             for sz=1:length(ph_num1)
                                 ph_num1(sz)=phage_aktiv(sz,randi(ph_num1(sz)));
                             end
                             ph_class=[ph_class; reshape(repmat(ph_num1,[1,phagenbinit])', [phagenbinit*length(ph_num1),1])];
                         end
                        
                        % update bacteria matrices:
                        positions = positions(idx,:);
                        dead_bacteria = dead_bacteria +bacnb-length(idx); 
                        bact_temp_state = bact_temp_state(idx);
                        bact_class=bact_class(idx);
                        bact_prop = bact_prop(idx,:);
                        bact_id=bact_id(idx);
                        bacnb = size(positions,1);
                        phage_array=phage_array(idx,:);
                        CRISPR_array=CRISPR_array(idx,:);
                        
                        simsteps
                    end
                    death_cnt = 0; % set counter back 
                end
                
                % phage decays after ph_decay step:
                ph_age=ph_age+1; % phage ageing
                ph_id=ph_age<ph_decay;
                ph_age=ph_age(ph_id);
                ph_class=ph_class(ph_id);
                ph_positions=ph_positions(ph_id,:);

% store the variables of bacteria:
popnum(simsteps,1)=size(positions,1); % bacterial population number
popnum(simsteps,2)=size(ph_positions,1); % phage popultion number
popnum(simsteps,3)=length(find(ph_freq_free)); % free phage diversity
popnum(simsteps,4)=length(find(ph_freq_crispr)); % crispr phage diversity
popnum(simsteps,5)=length(find(ph_freq_temp)); % temperate phage diversity
popnum(simsteps,6)=length(find(ph_freq_crispr+ph_freq_temp)); % phages in the metagenome (crispr+temperate)
popnum(simsteps,7)=uptake_crispr; % number of CRISPR spacer uptakes in each steps
popnum(simsteps,8)=uptake_temp; % number of temperate phage insertions in each steps

for i=1:SPECIES_NUM % bacteria numbers
    popnum(simsteps,8+i)=size(positions(bact_class==i,:),1);
end
for i=1:no_phage % phage numbers
    popnum(simsteps,8+SPECIES_NUM+i)=size(ph_positions(ph_class==i,:),1);
end   

% set bact the counters:
uptake_crispr=0;
uptake_temp=0;

% supplementary plotting step:
if mod(simsteps,community_development_step)==0
    histogrammainphage(phage_array,CRISPR_array, no_phage, no_crispr, simsteps);
    genomes(CRISPR_array, phage_array,simsteps);  
    Matrix(CRISPR_array,bact_class,phage_array,SPECIES_NUM,no_phage,simsteps);
end

end %end of the iterations
elapsedtime=toc(tstart)/60 %simulation time in minutes

%%
% visualisations:
if PLOT==0

    figure;
    hold on;
    findfigs;
    mainaxes = gca;
    set(mainaxes, 'ylim', [-8.5,4.5]);
    set(mainaxes, 'xlim', [-4.5,4.5]);
   colormap([linspace(0,1,no_phage+1)', linspace(1,0,no_phage+1)', 0*ones(no_phage+1,1)]);
   color_type = distinguishable_colors(no_phage);

       scatter(mainaxes,linspace(-4,4,scale_crispr), -8*ones(1,scale_crispr), (1:scale_crispr).^2*szorzo, linspace(scale_ph,0,scale_crispr), 'filled'); % plot CRISPR circles
       for i=1:no_phage % initial plot the phages
       plot(mainaxes,ph_positions(ph_class==i,1),ph_positions(ph_class==i,2), 'Color', color_type(i,:), 'MarkerFaceColor', color_type(i,:), 'marker', 'o','MarkerSize', 0.001, 'linestyle','none');
       end
       
       for i=1:SPECIES_NUM % change if diff number of phage and cell types  
       plot(mainaxes,positions(bact_class==i,1),positions(bact_class==i,2),'.', 'Color', color_type(i,:), 'MarkerFaceColor', color_type(i,:), 'marker', 'o','MarkerSize', 3,'linestyle','none'); % plot final bacteria
       end
      
       xlabel(['circle size is proportional to the number of CRISPR spacers (0-', num2str(scale_crispr-1), ')']);
       set(gca, 'XTick', [], 'YTick', []);
       hcb1=colorbar;
       title(hcb1, 'number of integrated temperate phages');
       print('Number of integrated temperate phages');
      
end
       
    % SUBPLOT FOR free phage frequencies:
    subaxes_phage= axes('Position', [.2, .19, .25, .07]);
    bar(subaxes_phage, (1:no_phage), ph_freq_free, 'b', 'EdgeColor', 'b');
    set(subaxes_phage, 'ylim', [0,1], 'xlim', [0,no_phage], 'Fontsize', 7);
    set(get(subaxes_phage,'XLabel'),'String','phages', 'Fontsize', 7);
    set(get(subaxes_phage,'YLabel'),'String','ferquency', 'Fontsize', 7);
    
    % SUBPLOT FOR temperate phage frequencies:
    subaxes_phagetemp= axes('Position', [.2, .31, .25, .07]);
    bar(subaxes_phagetemp, (1:no_phage), ph_freq_temp, 'r', 'EdgeColor', 'r');
    set(subaxes_phagetemp, 'ylim', [0,1], 'xlim', [0,no_phage], 'Fontsize', 7);
    set(get(subaxes_phagetemp,'XLabel'),'String','Temperate phages', 'Fontsize', 7);
    set(get(subaxes_phagetemp,'YLabel'),'String','Ferquency', 'Fontsize', 7);
    
    % SUBPLOT FOR CRISP spacer frequencies:
    subaxes_phagecrispr= axes('Position', [.5, .31, .25, .07]);
    bar(subaxes_phagecrispr, (1:no_phage), ph_freq_crispr, 'g', 'EdgeColor', 'g');
    set(subaxes_phagecrispr, 'ylim', [0,1], 'xlim', [0,no_phage], 'Fontsize', 7);
    set(get(subaxes_phagecrispr,'XLabel'),'String','CRISPR spacers', 'Fontsize', 7);
    set(get(subaxes_phagecrispr,'YLabel'),'String','ferquency', 'Fontsize', 7);
    
    print(['finalstate_',num2str(simsteps) ],'-djpeg');
   
mind2=nan(size(CRISPR_array));
for i=1:size(CRISPR_array, 1)
[van, hol]=ismember(CRISPR_array(i,:), phage_array(i,:));
mind2(i,hol~=0)=CRISPR_array(i,hol~=0);
end

% final state plots:
 figure
 subplot(1,2,1)
 hold on; 
 plot(popnum(:,1), 'LineStyle', '-', 'Color', [0.5 0.5 0.5], 'Display', 'All bacteria', 'LineWidth', 4);
 for i=1:SPECIES_NUM
 plot(popnum(:,8+i), 'Color', color_type(i,:), 'Display', ['Bacteria ', num2str(i)], 'LineWidth', 2);
 end
 xlabel('steps');
 ylabel('cell number');
 ylim([0, POPULATION_SIZE]);
 legend;
 subplot(1,2,2)
 hold on; 
 plot(popnum(:,2), 'LineStyle', '-', 'Color', [0.5 0.5 0.5], 'Display', 'All phage', 'LineWidth', 4);
 for i=1:no_phage
 plot(popnum(:,8+SPECIES_NUM+i), 'Color', color_type(i,:), 'Display', ['Phage ', num2str(i)], 'LineWidth', 2);
 end
 xlabel('steps');
 ylabel('phage number');
 legend;
 print(['bact_phage_num_',num2str(simsteps) ],'-djpeg');
 
  histogrammainphage(phage_array, CRISPR_array, no_phage, no_crispr, simsteps);
  genomes(CRISPR_array, phage_array,simsteps);   
  Matrix(CRISPR_array,bact_class,phage_array,SPECIES_NUM,no_phage,simsteps);

%%
save(bacteria_population_envf, 'popnum', 'free_freq', 'temp_freq', 'crispr_freq'); %save output data to .mat file 