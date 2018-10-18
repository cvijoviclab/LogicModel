function [PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, placeholders] = ...
   runUntilSteadyState(PromSiteIn, MetabolitesIn, MisclIn, Snf1pwIn, R2S3pwIn, PKApwIn, TORpwIn, placeholdersIn, activeCrosstalks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run model until it is in an attractor state given starting conditions
% specified in PromSiteIn, MetabolitesIn, MisclIn, Snf1pwIn, R2S3pwIn,
% PKApwIn, TORpwIn and placeholders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set start values

Snf1pw = Snf1pwIn;
Miscl = MisclIn;
PKApw = PKApwIn;
R2S3pw = R2S3pwIn;
TORpw = TORpwIn;
Metabolites = MetabolitesIn;
PromSite = PromSiteIn;
placeholders(1) = placeholdersIn(1);
placeholders(2) = placeholdersIn(2);

% RUN UNTIL ATTRACTOR STATE (= NOTHING CHANGES ANYMORE)

iteration = 1;

while true
    
    % save old values to compare later if something has changed
    Snf1pwOld = Snf1pw;
    MisclOld = Miscl;
    PKApwOld = PKApw;
    R2S3pwOld = R2S3pw;
    TORpwOld = TORpw;
    MetabolitesOld = Metabolites;
    PromSiteOld = PromSite;
    
    
    % glucose
    if (Metabolites{1,2} >= 1) && (Miscl{1,2} == 1 )
        Metabolites{1,3} = 1 ;
        Metabolites{2,2} = 1;
        Metabolites{2,3} = 1;
    end
    %import by Hexosetransporters PMID: 9299703
    
    if (Metabolites{1,2} == 0)
        Metabolites{1,3} = 0;
        Metabolites{2,2} = 0;
        Metabolites{2,3} = 0;
        Metabolites{3,2} = 0;
        Metabolites{3,3} = 0;
    end 
    
    %nitrogen
    if (Metabolites{4,2} == 1)
        Metabolites{4,3} = 1; % in cytosol
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Snf1 pathway %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%Reg1 dephosphorelation under glucose starvation by unknown protein
    %%%Xxx3
    if (Snf1pw{2,2} ==1) && (Miscl{6,2} == 1) && (Metabolites{1,2} == 0)
        Snf1pw{2,4} = 0;
    end
    
    %%%Glc7 dephosphorelation under glucose starvation by unknown protein
    %%%Xxx3
    if (Snf1pw{3,2} ==1) && (Miscl{6,2} == 1) && (Metabolites{1,2} == 0)
        Snf1pw{3,4} = 0;
    end
    
    %%%Snf1 phosphorylation/dephosphorylation
    if ((Snf1pw{1,2} ==1) && ((Snf1pw{4,2}==1)||(Snf1pw{5,2}==1)||(Snf1pw{6,2}==1)))
        Snf1pw{1,4} = 1 ;
    end %constituative phosphorylation of Snf1 by Sak1, Elm1 and Tos3 PMID:15831494,PMID:17991748

    if (Snf1pw{1,2} ==1) && (Snf1pw{2,2} ==1) && (Snf1pw{3,2} ==1) && (Snf1pw{2,4} == 1)
        Snf1pw{1,4} = 0; 
    end%dephosphorylation of Snf1 by Glc7 and phosphorylated Reg1 PMID:15831494

    %%%%Mig1 phosphorylation/dephosphorylation%%%%
    if (Snf1pw{1,2}==1) && (Snf1pw{10,2}==1) && ((Snf1pw{7,2}==1) || (Snf1pw{8,2}==1) || (Snf1pw{9,2}==1)) && (Snf1pw{1,4}==1)
        Snf1pw{11,4}=1;
    end%The SNF1 kinase complex beta and gamma subunits are needed for active Snf1 PMID:2557546  and PMID10990457
    
    if (Snf1pw{11,2} ==1) && (Snf1pw{2,2} ==1) && (Snf1pw{3,2} ==1) && (Snf1pw{2,4} == 1)
        Snf1pw{11,4}= 0;
    end% Mig1 dephosphorylation by Glc7 and phosphorylated Reg1 PMID: 28854669

    if (Snf1pw{11,2} == 1) && (Snf1pw{11,4} == 0)
        Snf1pw{11,3} = 2;
        Snf1pw{11,6} = 1;
    end%unphosphorylated Mig1 is nuclear and DNAbound PMID:9832517
    
    if (Snf1pw{11,2} == 1) && (Snf1pw{11,4} == 1)
        Snf1pw{11,3} = 1;
        Snf1pw{11,6} = 0;
    end%phosphorylated Mig1 is cytosolar 

    if (Snf1pw{11,6} == 1) && (Miscl{2,2}==1) && (Miscl{3,2}==1)
      Miscl{2,6}=1;
      Miscl{3,6}=1;
    end%Ssn6-Tup1 is recruited by Mig1 PMID:7724528

    if (Snf1pw{11,6} == 1) && (Miscl{2,6}==1) && (Miscl{3,6}==1)
      PromSite{1,2} = 0; %binding of Mig1 to SUC2 (repression)  PMID:9832517
      PromSite{2,2} = 0; %binding (repression) of Mig1 to GAL1 and GAL4 PMID: 9973625 
      PromSite{3,2} = 0;
    else
      PromSite{1,2} = 1;
      PromSite{2,2} = 1;
      PromSite{3,2} = 1;  
    end%Mig1 represses MAL61, MAL62, Mal63 PMID:8529272 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%Rgt2/Snf3 pathway

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%glucose sensing and Yck1/2 recruitment
    if (R2S3pw{3,2}==1) && ((R2S3pw{1,2}==1) && (Metabolites{1,2} == 2)) || ((R2S3pw{2,2}==1) && (Metabolites{1,2} == 1))
        R2S3pw{3,3}=0;
    end%Rgt2 is a glucose sensor PMID: 8901598, and acts through Yck1 PMID: 14755054
    if (R2S3pw{4,2}==1) && ((R2S3pw{1,2}==1) && (Metabolites{1,2} == 2)) || ((R2S3pw{2,2}==1) && (Metabolites{1,2} == 1))
        R2S3pw{4,3}=0;
    end%Rgt2 is a glucose sensor PMID: 8901598, and acts through Yck2 PMID: 14755054
    if (R2S3pw{3,2}==1) && ((R2S3pw{1,2}==1) && (Metabolites{1,2} >= 1)) || ((R2S3pw{2,2}==1) && (Metabolites{1,2} >= 1))
        R2S3pw{3,3}=0;
    end%Snf3 is a glucose sensor PMID: 8901598, and acts through Yck1 and Yck2PMID: 020594
    if (R2S3pw{4,2}==1) && ((R2S3pw{1,2}==1) && (Metabolites{1,2} >= 1)) || ((R2S3pw{2,2}==1) && (Metabolites{1,2} >= 1))
        R2S3pw{4,3}=0;
    end%Snf3 is a glucose sensor PMID: 8901598, and acts through Yck1 and Yck2PMID: 020594

    %%%YCK phosphorylation of Mth1 Std1
    %Yck1/2 is responsible for phosphorylation of Mth1 and Std1 PMID: 14755054
    if(R2S3pw{3,2}==1) && (R2S3pw{3,3}==0) && (R2S3pw{5,2}==1)
        R2S3pw{5,4}=1;
    end
    if(R2S3pw{4,2}==1) && (R2S3pw{4,3}==0) && (R2S3pw{5,2}==1)
        R2S3pw{5,4}=1;
    end
    if(R2S3pw{3,2}==1) && (R2S3pw{3,3}==0) && (R2S3pw{6,2}==1)
        R2S3pw{6,4}=1;
    end
    if(R2S3pw{4,2}==1) && (R2S3pw{4,3}==0) && (R2S3pw{6,2}==1)
        R2S3pw{6,4}=1;
    end%Mth1 and Std1 are activated through YCK PMID: 020594 and 14755054

    %%%Mth1 Std1 dephosphorylation of Rgt1
    if(R2S3pw{7,2}==1) && (R2S3pw{5,2}==1) && (R2S3pw{5,4}==0)
        R2S3pw{7,4} = 0;
    end
    if(R2S3pw{7,2}==1) && (R2S3pw{6,2}==1) && (R2S3pw{6,4}==0)
        R2S3pw{7,4} = 0;
    end%inactivation (phosphorylation) of Mth1 and Std1 leads to the hyperphosphorylation of Rgt1 Pubmed ID: 12925759
    
    %%%Mth1 dephosphorelation under glucose withdrawal by Psy2 and Pph3
    if (R2S3pw{5,2} ==1) && (R2S3pw{8,2} ==1) && (R2S3pw{9,2} ==1) && (Metabolites{1,2} == 0)
        R2S3pw{5,4}= 0;
    end% Mth1 dephosphorylation by Psy2 and Pph3 PMID: 24277933
    
    %%%Sdt1 dephosphorelation under glucose withdrawal by Xxx1
    if (R2S3pw{6,2} ==1) && (Miscl{4,2} ==1) && (Metabolites{1,2} == 0)
        R2S3pw{6,4}= 0;
    end% Sdt1 dephosphorylation by Xxx1
    
    %%%Rgt1 dephosphorelation under glucose withdrawal by Xxx1
    if (R2S3pw{7,2} ==1) && (Miscl{4,2} ==1) && (Metabolites{1,2} == 0)
        R2S3pw{7,4}= 0;
    end% Rgt1 dephosphorylation by Xxx1
    
    %%%repression by Rgt1
    if(R2S3pw{7,2} == 1) && (R2S3pw{7,4} == 0)
        R2S3pw{7,6} = 1;
    end% Dephosporylation of Rgt1 returns its DNA binding capacity PMID: 12861007

    if(R2S3pw{7,6} == 1) && (Miscl{2,2}==1)  && (R2S3pw{7,4} == 0) && (Miscl{3,2}==1)
      Miscl{2,6} = 1;
      Miscl{3,6} = 1;
    end

    if(R2S3pw{7,6} == 1) && (R2S3pw{7,4} == 0) && (Miscl{2,6}==1) && (Miscl{3,6}==1)
        PromSite{4,2} = 0;
        PromSite{7,2} = 0;
        PromSite{8,2} = 0;
    else
        PromSite{4,2} = 1;
        PromSite{7,2} = 1;
        PromSite{8,2} = 1;
    end
    %repression of Hxt1, Hxt2, Hxt3 and Hxt4 by Rgt1 under no and low glucose conditions PMID:7862149, 12527758, 15705057

    %%% Rgt1 as activator
    if((R2S3pw{7,2} == 1) && (R2S3pw{7,4} == 1))
        R2S3pw{7,6} = 1;
    end
    if((R2S3pw{7,6} == 1) && (R2S3pw{7,4} == 1))
        PromSite{7,2} = 1;
    end%Under glucose condition Rgt1 is hyperphosphorylated and acts as a activator PMID: 12527758
    %review on which hexostransporters are expressed onder which conditions     PMID26205245 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%AMPK pathway

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%Gpr1 branch
    if(PKApw{1,2} == 1) && (PKApw{2,2} == 1) && (Metabolites{1,2} >0)
        PKApw{2,5} = 2;
    end% Gpa2 is the beta subunit of Gpr1 and Gpa2 needs to be be associated with GTP to interact with Cyr1PMID:17983752

    if(PKApw{2,2} == 1) && (PKApw{9,2} == 1) && (PKApw{2,5} == 2)
        PKApw{9,3} = 0;
    end% Gpa2 is the beta subunit of Gpr1 and Gpa2 needs to be be associated with GTP to interact with Cyr1PMID:17983752

    %%%Ras branch
    if(Metabolites{1,3} == 1) && (PKApw{3,2} == 1) && ((PKApw{5,2} == 1) || (PKApw{6,2} == 1))
        PKApw{3,5} = 2;
    end%in presence of glucose Ras1 and Ras2 are associated with GTP by Cdc25 and Sdc25 PMID:9628870, 8206969

    if(PKApw{3,2} == 1) && (PKApw{3,5} == 2)
        PKApw{3,3} = 0;
    end%Ras2-GTP is anchered to the cytoplasm side of the cell membrane PMID:8430318 

    if(PKApw{3,2} == 1) && (PKApw{9,2} == 1) && (PKApw{3,5} == 2)
        PKApw{9,3} = 0 ;
    end%Cyr1 localizes to the plasma membrane PMID:1875942

    if(Metabolites{1,3} == 1) && (PKApw{4,2} == 1) && ((PKApw{5,2} == 1) || (PKApw{6,2} == 1))
       PKApw{4,5} = 2;
    end%in presence of glucose Ras1 and Ras2 are associated with GTP by Cdc25 and Sdc25 PMID:9628870, 8206969

    if(PKApw{4,2} == 1) && (PKApw{4,5} == 2)
        PKApw{4,3} = 0; 
    end%Ras2-GTP is anchered to the cytoplasm side of the cell membrane PMID:8430318

    if(PKApw{4,2} == 1) && (PKApw{9,2} == 1) && (PKApw{4,5} == 2)
        PKApw{9,3} = 0; 
    end%Cyr1 localizes to the plasma membrane PMID:1875942

    if(Metabolites{1,3} == 0) && (PKApw{3,2} == 1) && ((PKApw{7,2} == 1) || (PKApw{8,2} == 1))
        PKApw{3,5} = 1; 
    end%Ras proteins are GDP-loaded by Ira1 and Ira2 PMID:15339905, 1668647 

    if(Metabolites{1,3} == 0) && (PKApw{2,2} == 1) && ((PKApw{7,2} == 1) || (PKApw{8,2} == 1))
        PKApw{4,5} = 1;
    end%Ras proteins are GDP-loaded by Ira1 and Ira2 PMID:15339905, 1668647

    %%% convert on Adenylate cyclas
    if(PKApw{9,2} == 1) && (PKApw{9,3} == 0) && (Metabolites{2,2} == 1)
        Metabolites{3,2} = 1;
        Metabolites{3,3} = 1;
    end%Glucose addition to yeast causes a spike in cAMP PMID:9628870

    %%%PKA part
    if(PKApw{13,2} ==1) && ((PKApw{10,2} ==1)||(PKApw{11,2} ==1)||(PKApw{12,2}==1)) && (Metabolites{3,2} == 1)
        placeholders(1) = 1;
    else
        placeholders(1) = 0;
    end%placeholder for cAMP binding to the PKA complex

    %%%Downstream PKA
    if(PKApw{14,2}==1 && placeholders(1) == 1)
        PKApw{14,4} = 1;
    end%active PKA inactivates Rim15 by phosphorylation PMID: 16759348,15661010
    
    if(PKApw{15,2}==1 && placeholders(1) == 1)
        PKApw{15,4} = 1;
    end%active PKA inactivates Yak1 by phosphorylation PMID: 21255108

    if(PKApw{14,2}==1 && (Metabolites{1,2} == 0))
        PKApw{14,4} = 0;
    end% under glucose starvation Rim15 gets dephosphorylated by Xxx2
    
    if(PKApw{15,2}==1 && (Metabolites{1,2} == 0))
        PKApw{15,4} = 0;
    end% under glucose starvation Yak1 gets dephosphorylated by Xxx2
    
    if(PKApw{14,2}==1 && PKApw{16,2}==1 && PKApw{14,4}==0)
        PKApw{16,4}=1;
    end%dephoshporylated Rim15 phosphorylates Msn2 PMID: 24140345

    if(PKApw{14,2}==1 && PKApw{17,2}==1 && PKApw{14,4}==0)
        PKApw{17,4}=1;
    end%dephoshporylated Rim15 stimulates Msn4 PMID: 24140345
    
    if(PKApw{15,2}==1 && PKApw{16,2}==1 && PKApw{15,4}==0)
        PKApw{16,4}=1;
    end%Yak1 phosphorylates Msn2 PMID: 18793336 

    if(PKApw{15,2}==1 && PKApw{17,2}==1 && PKApw{15,4}==0)
        PKApw{17,4}=1;
    end%Yak1 phosphorylates Msn4 PMID:18793336 
    
    if(placeholders(1) == 1 && PKApw{16,2}==1)
        PKApw{16,3}=1;
        PKApw{16,6}=0;
    end%PKA represses nuclear localization of Msn2 PMID: 9472026, 16281053
    
    if(placeholders(1) == 0 && PKApw{16,2}==1)
        PKApw{16,3}=2;
    end%PKA represses nuclear localization of Msn2 PMID: 9472026, 16281053

    if(placeholders(1) == 1 && PKApw{17,2}==1)
        PKApw{17,3}=1;
        PKApw{17,6}=0;
    end%PKA represses nuclear localization of Msn4 PMID: 9472026
    
    if(placeholders(1) == 0 && PKApw{17,2}==1)
        PKApw{17,3}=2;
    end%PKA represses nuclear localization of Msn4 PMID: 9472026
    
    if PKApw{14,2}==1 && PKApw{14,4}==1 && PKApw{18,2}==1
        PKApw{18,4} = 0;
        PKApw{18,6} = 0;
    end% inactive Rim15 deactivates Gis1 (connected to Igo1/2, Cdc55/Pph21/Tpd3 complex) PMID: 23273919
    
    if PKApw{16,2}==1 && Miscl{7,2} == 1 && Metabolites{1,2} >= 1
        PKApw{16,4} = 0;
    end% unknown protein Xxx4 dephoshorelates Msn2 in glucose
    
    if PKApw{17,2}==1 && Miscl{7,2} == 1 && Metabolites{1,2} >= 1
        PKApw{17,4} = 0;
    end% unknown protein Xxx4 dephoshorelates Msn4 in glucose
      

    if(PKApw{14,2}==1 && PKApw{18,2}==1 && PKApw{14,4}==0)
        PKApw{18,4}=1;
    end%Rim15 stimulates Gis1 PMID: 15300954,10835355 

    if(PKApw{16,2}==1 && PKApw{16,3}==2 && PKApw{16,4}==1)
        PKApw{16,6}=1;
        PromSite{5,2} = 1;
    else
        PromSite{5,2} = 0;
    end%Msn2 activates STRE genes PMID:8650168, 8641288

    if(PKApw{17,2}==1 && PKApw{17,3}==2 && PKApw{17,4}==1)
        PKApw{17,6}=1;
        PromSite{5,2} = 1;
    else
        PromSite{5,2} = 0;
    end%Msn4 activates STRE genes PMID:8650168, 8641288

    if(PKApw{18,2}==1 && PKApw{18,4}==1 && PKApw{18,3}==2 )
        PKApw{18,6}=1; 
        PromSite{6,2} = 1; 
    else
        PromSite{6,2} = 0;
    end%Gis1 activates PDS genes PMID:22363679
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%TOR pathway

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% in nitrogen starvation and glucose Par32 gets phosphorylated by
    %%% Gtr1 and Npr2, PMID:???25085507
    if (Metabolites{4,2} == 0) && (Metabolites{1,2} == 1) && (TORpw{13,2}== 1) && (TORpw{14,2}== 1) && (TORpw{15,2}== 1)
        TORpw{13,4} = 1;
    end
    
    %%% TORC1 gets activated through AA and Eco complex if nitrogen is
    %%% available and deactivated if not
    if (Metabolites{4,2} == 1) && ((TORpw{1,2}== 1)||(TORpw{2,2}== 1)) && ...
            ((TORpw{3,2}== 1) && (TORpw{3,4}== 0)) && (TORpw{4,2}== 1) && (TORpw{5,2}== 1)
        placeholders(2) = 1; %%% active TORC1
    else
        placeholders(2) = 0; %%% inactive TORC1
    end

    %%% Sch9 phosphorylation by active TORC1, PMID: 22964838
    if (placeholders(2) == 1) && (TORpw{12,2} ==1)
        TORpw{12,4} = 1 ;
    end 
    
    %%% Tap42 phosphorylation by active TORC1, PMID: 22964838, 22174183
    if (placeholders(2) == 1) && (TORpw{6,2} ==1)
        TORpw{6,4} = 1;
    end

    %%% then Tap42 binds to catalytic subunits of PP2A
    %%% (Php21/22/3 and Sit4), is bound to the vacuolar membrane,
    %%% and can therefore not dephosphorylate Gln3, PMID: 12820961, 25085507  
    if ((TORpw{6,2} ==1) && (TORpw{6,4} == 1)) && (TORpw{7,2} ==1) && (TORpw{8,2} ==1) 
        if TORpw{10,2} == 1
            TORpw{10,4} = 1;
        end
        if (TORpw{11,2} == 1)
            TORpw{11,3} = 1;
            TORpw{11,6} = 0;
            TORpw{11,4} = 1;
        end
    end
    
    %%% Sch9 dephosphorylated by deactive TORC1, PMID: 22964838
    if (placeholders(2) == 0) && (TORpw{12,2} ==1)
        TORpw{12,4} = 0;
    end 

    %%% Tap42 dephosphorylated by PP2A and no TORC1, PMID: 22964838,
    %%% 22174183,???24738657
    if (placeholders(2) == 0) && (TORpw{6,2} ==1) && (TORpw{7,2} ==1)
        TORpw{6,4} = 0;
    end
    
    %%% then Tap42 is released to cytosol, the complex with ...
    %%% PP2A(Php21/22/3 and Sit4) is not present, PMID: 12820961
    %%% if there is phosphorylated Par32 that
    %%% activates the PP2A branch, PMID: 20093466
    
    %%% Tip41 blocks Tap42 to bind to Sit4/PP2A, PMID: 11741537
    %%% free Sit4/PP2A promotes dephosphorylation of Gln3, Ure2
    %%% Gln3 and Ure2 release their binding
    
    %%% Ure2 binds to both hyperphosphorylated, phosphorylated and unphosphorylated Gln3, but
    %%% only if it phosphorylated itself, otherwise it will be in a complex
    %%% with Mks1, and can therefore not bind to Gln3, PMID: 10940301
    
    %%% dephosphorylated Gln3 by TOR signaling moves to the nucleus, PMID: 22174183,22964838 
    
    if (TORpw{13,2} == 1) && (TORpw{13,4} == 1) && (TORpw{6,2} == 1) && ...
            (TORpw{6,4} == 0) && (TORpw{7,2} == 1) && (TORpw{8,2} == 1) && (TORpw{9,2} == 1)
        if TORpw{10,2} == 1
            TORpw{10,4} = 0;
        end
        if TORpw{11,2} == 1
            TORpw{11,4} = 0;
            TORpw{11,3} = 2;
            TORpw{11,6} = 1;
        end
    else %Gln3 is in the cytosol if the PP2A branch is inactivated and there is no Snf1 that would put it in again (no glucose, crosstalk(10) is on), has to be here: looping problem otherwise
        if TORpw{10,2} == 1
            TORpw{10,4} = 1;
        end
        if TORpw{11,2} == 1 && TORpw{10,2} == 1 && TORpw{10,4} == 1 && ~(Metabolites{1,2} == 0 && activeCrosstalks(10) == 1)
            TORpw{11,4} = 1;
            TORpw{11,3} = 1;
            TORpw{11,6} = 0;
        end
    end

    %%% cytosolic Gln3 cannot express NCR gene
    if (TORpw{11,2} == 1) && (TORpw{11,3} == 1)
        PromSite{9,2} = 0;
    else
        PromSite{9,2} = 1;
    end

    %%% nuclear Gln3 expresses NCR gene, PMID: 19104072, 20378536
    if (TORpw{11,2} == 1) && (TORpw{11,3} == 2) 
        PromSite{9,2} = 1;
    else
        PromSite{9,2} = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%Crosstalk

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, placeholders] = ...
    crosstalks(Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, placeholders, activeCrosstalks);
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%Stop criterium

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if nothing changes break the loop
    if isequal(Snf1pwOld, Snf1pw) && isequal(MisclOld, Miscl) && isequal(PKApwOld, PKApw) && ...
            isequal(R2S3pwOld, R2S3pw) && isequal(MetabolitesOld, Metabolites) && ...
            isequal(PromSiteOld, PromSite) && isequal(TORpwOld, TORpw)
        disp('Steady state reached');
        break;
    else
        disp(['Iteration ', num2str(iteration), ' completed']);
    end
    
    % if steady state is not reached after 100 iterations stop
    if iteration == 100
        disp('Something went wrong. Steady state is not reached after 100 iterations');
        disp(' ');
        break;
    end
    
    iteration = iteration + 1;
    
end

end