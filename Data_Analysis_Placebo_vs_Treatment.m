%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                         %        Data Analysis       %
                         %    Placebo vs Treatment    %
                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all

% Need to be in the files data
files_1 = dir('repro_d1*');
files_2 = dir('repro_d2*');

for i =1 :size(files_1,1)
    
    %% LOAD FILES
    tempfile_1(i) = load (files_1(i).name);
    tempfile_2(i) = load (files_2(i).name);

    
    %% REGRESSION FOR EACH CONDITION
    
R = regress(zscore(tempfile_1(i).data.reproduction(:,2)), [ones(120,1) zscore(tempfile_1(i).data.reproduction(:,1))]);
R_1 (i) = R(2,1);
R_t  = regress(zscore(tempfile_2(i).data.reproduction(:,2)), [ones(120,1) zscore(tempfile_2(i).data.reproduction(:,1))]);
R_2 (i) = R_t(2,1);


  if strcmp(tempfile_1(i).data.condition, 'placebo') 
      
       
       R_placebo (i)  = R_1 (i);
       R_panortaxin (i)  = R_2 (i);
  else 
   
        R_panortaxin (i)  = R_1(i);
        R_placebo (i)  = R_2(i);
  end
  
  
   %% EXTRACTION DATA
   
   % Extraction age 
   age(i) = tempfile_1(i).data.age;
   
   % Extraction gender
   % Male = 0
   % Female = 1
   
   if strcmp(tempfile_1(i).data.gender, 'male')
       
       gender (i) = 0;
       
   else 
       
       gender (i) = 1;
       
   end
end   
   

   % Extraction mean RT for each trials

for trial=400:200:1400
    for i=1:size(files_1,1)

     h = (trial/200)-1;
     
          if strcmp(tempfile_1(i).data.condition, 'placebo')

                M_placebo (h,i)  = mean(tempfile_1(i).data.reproduction(tempfile_1(i).data.reproduction(:,1)==trial, 2));
                M_panortaxin (h,i)  = mean(tempfile_2(i).data.reproduction(tempfile_2(i).data.reproduction(:,1)==trial, 2));

          else 

                M_panortaxin (h,i)  = mean(tempfile_1(i).data.reproduction(tempfile_1(i).data.reproduction(:,1)==trial, 2));
                M_placebo (h,i)  = mean(tempfile_2(i).data.reproduction(tempfile_2(i).data.reproduction(:,1)==trial, 2));

  
           end
     end 
 end
 
 
 
     M_placebo = [M_placebo]';
     M_panortaxin = [M_panortaxin]';
     
 % Mean total 
 
 for h=1:size(M_placebo,1)
     
     Mean_placebo_total (h) = mean(M_placebo(h,:));
     Mean_panortaxin_total  (h)= mean(M_panortaxin(h,:));
 end 
 
 

    % Summary of participant ADD GENDER
    
    % First column ==> age
    % Second column ==> gender (male = 0 // female = 1) 
    % Third column ==> ßeta condition placebo
    % Fourth column ==> mean RT placebo 
    % Fifth - tenth columns ==> mean RT for 400 - 1400 ms placebo
    % eleventh column ==> ßeta condition panortaxin
    % twelfth column ==> mean RT panortaxin 
    % thirteenth - eighteenth ==> mean RT for 400 - 1400 ms panortaxin
    
    
   summary = [[age]' [gender]' [R_placebo]' [Mean_placebo_total]' [M_placebo] [R_panortaxin]' [Mean_panortaxin_total]' [M_panortaxin] ];
   
   
   %% T-TEST ==> GENDER
   
   [h1,p1,ci1,t_test_mean_placebo_gender] = ttest2(summary(summary(:,2)==0,4), summary(summary(:,2)==1,4));
   
   %% PEARSON CORRELATION ==> AGE ~(MEAN PLACEBO)
   
   [rho, pval] =  corrcoef(summary(:,1),summary(:,4), 'alpha', 0.05);
   df = (size (summary(:,4),1))-2;
   
   
   %% T-TEST ==> ==> MEAN ~ CONDITION
   
   [h2,p2,ci2,t_test_mean_condition] = ttest(summary(:,4), summary(:,12));
   
    %% T-TEST ==> ==> ßETA ~ CONDITION
   
   [h3,p3,ci,t_test_beta_condition] = ttest(summary(:,3), summary(:,11));
   
   %% DEMOGRAPHIC INFORMATION 
   % Male 
   age_mean_male = mean(summary(summary(:,2)==0, 1));
   age_std_male = std(summary(summary(:,2)==0, 1));
   total_male = sum(summary(:,2)==0);
   
   % Female 
   age_mean_female = mean(summary(summary(:,2)==1, 1));
   age_std_female = std(summary(summary(:,2)==1, 1));
   total_female = sum(summary(:,2)==1);
   
   % Total  
   age_mean_total = mean(summary(:,1));
   age_std_total = std(summary(:,1));
   age_total = sum(total_male)+sum(total_female);
   
   %% STATISTICAL INFORMATION
   % Mean RT condition placebo ~ Gender
   mean_condition_placebo_male = mean(summary(summary(:,2)==0,4));
   mean_condition_placebo_female = mean(summary(summary(:,2)==1,4));
   mean_condition_placebo = [mean_condition_placebo_male mean_condition_placebo_female];
   std_condition_placebo_male = std(summary(summary(:,2)==0,4));
   std_condition_placebo_female = std(summary(summary(:,2)==1,4));
   std_condition_placebo = [std_condition_placebo_male std_condition_placebo_female];
 
   % Mean RT ~ Condition
   mean_RT_placebo = mean(Mean_placebo_total);
   mean_RT_panortaxin = mean(Mean_panortaxin_total);
   mean_RT = [mean_RT_placebo mean_RT_panortaxin];
   std_RT_placebo = std(Mean_placebo_total);
   std_RT_panortaxin = std(Mean_panortaxin_total);
   std_RT = [std_RT_placebo std_RT_panortaxin];
   
    % ßeta RT ~ Condition
    beta_placebo = mean(summary(:,3));
    beta_panortaxin = mean(summary(:,11));
    beta_condition = [beta_placebo beta_panortaxin];
    std_beta_placebo = std(summary(:,3));
    std_beta_panortaxin = std(summary(:,11));
    std_beta_condition = [std_beta_placebo std_beta_panortaxin];
    
 %% PRINT 
 
 
   fprintf('Demographic Statistics\n\n');

fprintf('Age male:       Mean = %2f          Standard deviation = %2f          Sample size = %d\n\n',age_mean_male,age_std_male,total_male);
fprintf('Age female:     Mean = %2f          Standard deviation = %2f          Sample size = %d\n\n',age_mean_female,age_std_female,total_female);
fprintf('Age total:      Mean = %2f          Standard deviation = %2f          Sample size = %d\n\n',age_mean_total,age_std_total,age_total);

fprintf('                                \n\n ');
fprintf('                                \n\n ');
   fprintf('Inferential Statistics\n\n');
   
fprintf('Pearson Correlation : Age - Mean RT Pearson  \n\n');
fprintf('            Beta = %2f                     p-value =    %2f                \n\n',rho(1,2), pval(1,2));
fprintf('              df = %2f                          \n\n',df);
fprintf('                                \n\n ');

fprintf('T-test :  Mean RT  ~ Gender  \n\n');
fprintf('          t-stat = %2f                     p-value =  %2f                \n\n',t_test_mean_placebo_gender.tstat, p1);
fprintf('             df  = %2f                      Std   =  %2f            \n\n',t_test_mean_placebo_gender.df,t_test_mean_placebo_gender.sd);
fprintf('      Mean male  = %2f                 Std male   =  %2f            \n\n',mean_condition_placebo_male,std_condition_placebo_male);
fprintf('    Mean female  = %2f                Std female  =  %2f            \n\n',mean_condition_placebo_female,std_condition_placebo_female);

fprintf('                                \n\n ');
fprintf('Paired t-test :  Mean RT  ~ Condition  \n\n');
fprintf('           t-stat = %2f                    p-value =  %2f                \n\n',t_test_mean_condition.tstat, p2);
fprintf('              df  = %2f                     Std   =  %2f            \n\n',t_test_mean_condition.df,t_test_mean_condition.sd);
fprintf('    Mean Placebo  = %2f             Std Placebo   =  %2f            \n\n',mean_RT_placebo,std_RT_placebo);
fprintf(' Mean Panortaxin  = %2f           Std Panortaxin  =  %2f            \n\n',mean_RT_panortaxin,std_RT_panortaxin);


fprintf('                                \n\n ');
fprintf('Paired t-test :  ßeta RT  ~ Condition  \n\n');
fprintf('           t-stat = %2f                    p-value =   %2f                \n\n',t_test_beta_condition.tstat, p3);
fprintf('              df  = %2f                     Std   =   %2f            \n\n',t_test_beta_condition.df,t_test_beta_condition.sd);
fprintf('    ßeta Placebo  = %2f               Std Placebo   =  %2f            \n\n',beta_placebo,std_beta_placebo);
fprintf(' ßeta Panortaxin  = %2f             Std Panortaxin  =  %2f            \n\n',beta_panortaxin,std_beta_panortaxin);



   %% PLOT 
   %% Plot mean ==> placebo ~ gender
   subplot (3,2,1)
   se_condition_placebo_male = std_condition_placebo_male/sqrt(total_male);
   se_condition_placebo_female = std_condition_placebo_female/sqrt(total_female);
   se_condition_placebo = [se_condition_placebo_male se_condition_placebo_female];
   errorbar ( mean_condition_placebo, se_condition_placebo, '.');
   str_1 = {'Male', 'Female'};
   hold on
   set(gca, 'XTickLabel',str_1, 'XTick',1:numel(str_1));
   ylim([898 901]);
   xlim([0.5 2.5]);
   ylabel('Time (ms)')
   title('A')


 
   %% Scatter plot ==> mean RT placebo ~ age 

     subplot (3,2,2)
     scatter (summary(:,1), summary(:,4), 3,'k', '+');
     xlabel('Age')
     xlim([18 55])
     ylim([870 930])
     ylabel('Time (ms)')
     title('B')


   %% Plot mean ==> mean RT ~ condition
   
   subplot (3,2,3)
   mean_RT = [mean(Mean_placebo_total) mean(Mean_panortaxin_total)];
   se_RT = [std_RT/sqrt(300)];
   str_1 = {'Placebo', 'Panortaxin'};
   hold on
   
   errorbar(mean_RT,se_RT, '.');
   set(gca, 'XTickLabel',str_1, 'XTick',1:numel(str_1));
   ylim([865 905]);
   xlim([0.5 2.5]);
   ylabel('Time (ms)')
   title('C')
   
   %% Plot Beta ~ condition
   
    subplot (3,2,4)   
    
    beta_condition = [mean(summary(:,3))  mean(summary(:,11))];
    se_beta_condition = [std_beta_placebo/sqrt(300) std_beta_panortaxin/sqrt(300)];
    str_3 = {'Placebo', 'Panortaxin'};
    hold on
    errorbar(beta_condition,se_beta_condition,'.')
    set(gca, 'XTickLabel',str_3, 'XTick',1:numel(str_3));
    ylim([0.98 1]);
    xlim([0.5 2.5]);
    ylabel('Correlation coefficient')
    title('D')
    
    

   %% Plot mean ==> mean RT for each condition ~ trials
   % Mean and std of ßeta for each condition and trial
    for i = 1:6
        mean_RT_trials (i,1:2) = [mean(summary(:,(4+i))) mean(summary(:,(12+i)))];
        std_RT_trials (i,1:2) = [std(summary(:,(4+i))) std(summary(:,(12+i)))];
        
    end
    
  subplot (3,2,[5 6]) 
  errorbar(mean_RT_trials, (std_RT_trials/sqrt(300)))
  str_trials = {'400', '600', '800', '1000', '1200', '1400'};
  set(gca, 'XTickLabel',str_trials, 'XTick',1:numel(str_trials));
  legend({'Placebo';'Panortaxin'},'Location', 'northwest')
  xlim([0.5 6.5])
  ylim([0 1500])
  ylabel('Time (ms)')
  xlabel('Trials')
  title('E')

 
    saveas(gcf,'myfigure.pdf')
 
    
    
   
 

    
       
   


 
