function [Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, placeholders] = ...
    crosstalks(Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, placeholders, activeCrosstalk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, placeholders,
% is the data of the current state as input and the updated one as output

% activcrosstalk is a array with length n depending on how many crosstalks
% are added. Each position in the array is either 1 or 0, depending if the
% crosstalk is active or not

% this part of the model is here to simulate the crosstalk we have found in
% the literature.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Crosstalk 1%%%
    if ((activeCrosstalk(1)==1) && (Snf1pw{1,2} ==1) && (Snf1pw{1,4} == 0) && (R2S3pw{5,2}==1))
        R2S3pw{5,4}=1;
    end
    if ((activeCrosstalk(1)==1) && (Snf1pw{1,2} ==1) && (Snf1pw{1,4} == 0) && (R2S3pw{6,2}==1))
        R2S3pw{6,4}=1;
    end%Inactive Snf1 prevent degradation of Mth1 and Std1. PubMedID: 17586499, 16361229
    
%%%Crosstalk 2%%%
    if ((activeCrosstalk(2)==1) && (Snf1pw{1,2}==1) && (R2S3pw{5,2}==1) && (R2S3pw{5,4}==0))
        Snf1pw{1,4}=1;

    end%Std1 stimulates the Snf1 kinase activity. PubMedID: 8114728,12618390, 12220226
    
%%%Crosstalk 3%%%   
    if ((activeCrosstalk(3)==1) && (R2S3pw{3,2}==1) && (Snf1pw{2,2} ==1) && (Snf1pw{3,2} ==1) && (Metabolites{1,2} >= 1))
        R2S3pw{3,3}=0;
    end %Reg1-Glc7 acts as an upstream activator of Yck1. PubMedID: 16361229
    if ((activeCrosstalk(3)==1) && (R2S3pw{4,2}==1) && (Snf1pw{2,2} ==1) && (Snf1pw{3,2} ==1) && (Metabolites{1,2} >= 1))
        R2S3pw{4,3}=0;
    end %Reg1-Glc7 acts as an upstream activator of Yck2. PubMedID: 16361229 

%%%Crosstalk 4%%%      
    if (activeCrosstalk(4)==1) && (Snf1pw{4,2}==1) && (PKApw{13,2} ==1) && ((PKApw{10,2} ==1)||(PKApw{11,2} ==1)||(PKApw{12,2}==1)) && (placeholders(1) == 1)
       Snf1pw{4,4}=1;
    end %PKA complex phosphorylates Sak1. PubMedID:22140226

%%%Crosstalk 5%%%   
    if (activeCrosstalk(5)==1) && (Snf1pw{1,2} ==1) && (PKApw{13,2} ==1) && ((PKApw{10,2} ==1)||(PKApw{11,2} ==1)||(PKApw{12,2}==1)) && (placeholders(1) == 1)
        Snf1pw{1,4} = 0; 
    end %PKA complex negatively regulates the Snf1 pathway (Sak1 independent). PubMedID:22140226

%%%Crosstalk 6%%% 
    if (activeCrosstalk(6)==1) && (Snf1pw{1,2} ==1) && (Snf1pw{1,4} ==1)
        PKApw{16,4}=1;
    end %Snf1 can phosphorylate Msn2. PubMedID: 16281053 

%%%Crosstalk 7%%% 
    if (activeCrosstalk(7)==1) && (Snf1pw{2,2} ==1) && (Snf1pw{3,2} ==1) && ((PKApw{10,2} ==1)||(PKApw{11,2} ==1)||(PKApw{12,2}==1)) && (placeholders(1) == 1)
        Snf1pw{2,4}=1;
        Snf1pw{3,4}=1;
    end %glucose activation of the cAMP-PKA (protein kinase A) pathway is required for glucose activation of PP1. PubMedID: 22290422

%%%Crosstalk 8%%% 
    if (activeCrosstalk(8)==1) && (Snf1pw{1,2} ==1) && (Snf1pw{1,4} ==1) && (PKApw{9,2} == 1)
        PKApw{9,4} = 1;
    end
    
    if (PKApw{9,4} == 1)
        PKApw{9,3} = 1;
    end
    %Snf1 deactivates Cyr1 by phosphorylation. PubMedID: 26309257

%%%Crosstalk 9%%%
    if((activeCrosstalk(9)==1) && (R2S3pw{7,2}==1) && (PKApw{13,2} ==1) && ((PKApw{10,2} ==1)||(PKApw{11,2} ==1)||(PKApw{12,2}==1)) && (placeholders(1) == 1))
        R2S3pw{7,4}=1;
    end%Bcy1 phosphorylates Rgt1 under high glucose conditions. PubMedID: 21748783, 23468525, 16844691

    
%%%Crosstalk 10%%%
    if (activeCrosstalk(10)==1) && (Snf1pw{1,2}==1) && (Snf1pw{1,4}==1) && ...
            (Metabolites{4,2} == 1) && (TORpw{11,2}==1)
        TORpw{11,4} = 1;
        TORpw{11,3} = 2;
        TORpw{11,6} = 1;
    end % under glucose starvation (but not nitrogen starvation) phos. Snf1 can phosphorelate Gln3 and locate it in the nucleus, it's a different phosphorelation site, and has no connection to the TOR mediated phosphorelation of Gln3, PMID: 11809814,15911613  
    
    if (activeCrosstalk(10)==1) && (Snf1pw{1,2}==1) && (Snf1pw{1,4}==1) && ...
            (TORpw{3,2}==1)
        TORpw{3,4} = 1;
        placeholders(2) = 0;
    elseif (activeCrosstalk(10)==1) && (TORpw{3,2}==1) && (Miscl{9,2}==1) && Metabolites{1,2} == 1
        TORpw{3,4} = 0;
    end % under glucose starvation phos. Snf1 phosphorelates Kog1, PMID: 28096180, otherwise it is phosphorelated by an unknown mechanism Xxx6
    
    if (activeCrosstalk(10)==1) && (Snf1pw{1,2}==1) && (Snf1pw{1,4}==1) && ...
            (TORpw{12,2}==1)
        TORpw{12,4} = 0;
    end % under glucose starvation phos. Snf1 (partly) dephosphorelates Sch9, PMID: 25085507
    
    if (activeCrosstalk(10)==1) && ((Snf1pw{1,2}==1) && (Snf1pw{1,4}==1)) && ...
            (TORpw{13,2}==1)
        TORpw{13,4} = 0;
    end % under glucose starvation phos. Snf1 dephosphorelates Par32 and so deactivates the PP2A branch, PMID: 25085507

%%%Crosstalk 11%%%
    if (activeCrosstalk(11)==1) && (Miscl{8,2} == 1) && (Metabolites{1,2} == 0) && ...
            (Metabolites{4,2} == 0) && (TORpw{13,2}==1)
        TORpw{13,4} = 0;
    end % under glucose and nitrogen starvation Xxx5 can dephosphorelate Par32 and so deactivate the PP2A branch, PMID: 25085507
    
%%%Crosstalk 12%%%
    if (activeCrosstalk(12)==1) && (TORpw{12,2}==1) && (TORpw{12,4}==1) && ...
            (PKApw{14,2}==1)
        PKApw{14,3} = 1;      
    end% active Sch9 represses Rim15 nuclear localization, PMID: 20395504
    
    if (activeCrosstalk(12)==1) && (placeholders(2) == 1) && ...
            (PKApw{14,2}==1)
        PKApw{14,3} = 1;      
    end% active Tor represses Rim15 nuclear localization by dephosphorylation, PMID: 20395504, 16759348, 14690612
       
    if (activeCrosstalk(12)==1) && (placeholders(2) == 1) && ...
            (PKApw{16,2}==1)
        PKApw{16,3} = 1; 
        PKApw{16,6} = 0;   
    end% active Tor promotes nuclear export of Msn2, PMID: 12093809
    
%%%Crosstalk 13%%%   
    if (activeCrosstalk(13)==1) && (TORpw{12,2}==1) && (TORpw{12,4}==1) && ...
            (PKApw{17,2}==1)
        PKApw{18,4} = 1;  
        PKApw{18,6} = 1;
    end% active Sch9 phosphorylates Gis1 independent of Rim15 PMID: 20395504


end