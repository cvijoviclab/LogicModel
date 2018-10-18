function [PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, ...
    placeholders] = initalizeModel()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of the model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%Target promotors
Name = [string('SUC2'),string('GAL'),string('MAL'),string('HXT2/3/4'),string('STRE'),string('PDS'),string('HXT1'),string('HXK2'), string('NCR')]';
active = ones(length(Name),1);
PromSite = table(Name , active);
PromSite{5,2} = 0;
PromSite{6,2} = 0;


%%%%Metabolites
Name = [string('Glc'),string('ATP'),string('cAMP'),string('NH3')]';
presence = ones(length(Name),1);
localization = zeros(length(Name),1);
Metabolites = table(Name, presence, localization);
Metabolites{2,3} = 1;
Metabolites{3,2} = 0;

%%%%Proteins shared through pathways
Name = [string('HXTs'),string('Ssn6'),string('Tup1'), string('Xxx1'), ...
    string('Xxx2'), string('Xxx3'), string('Xxx4'), string('Xxx5'), string('Xxx6')]';
presence = ones(length(Name),1);
localization = ones(length(Name),1);
phosphorylation = zeros(length(Name),1);
GXP = zeros(length(Name),1);
DNA_binding = zeros(length(Name),1);
Miscl = table( Name, presence , localization, phosphorylation, GXP, DNA_binding);
Miscl{1,3} = 0;
Miscl{2,3} = 2;
Miscl{3,3} = 2;

%%%%Snf1 pathway
Name = [string('Snf1'),string('Reg1'),string('Glc7'),string('Sak1'),string('Tos3'),string('Elm1'),string('Sip1'),string('Sip2'),string('Gal83'),string('Snf4'),string('Mig1')]';
presence = ones(length(Name),1);
localization = ones(length(Name),1);
phosphorylation = zeros(length(Name),1);
GXP = zeros(length(Name),1);
DNA_binding = zeros(length(Name),1);
Snf1pw = table(Name , presence, localization, phosphorylation, GXP, DNA_binding);

%%%%Rgt2/Snf3 pathway
Name = [string('Rgt2'),string('Snf3'),string('Yck1'),string('Yck2'),string('Mth1'),...
    string('Std1'),string('Rgt1'),string('Pph3'),string('Psy2')]';
presence = ones(length(Name),1);
localization = ones(length(Name),1);
phosphorylation = zeros(length(Name),1);
GXP = zeros(length(Name),1);
DNA_binding = zeros(length(Name),1);
R2S3pw = table(Name , presence, localization, phosphorylation, GXP, DNA_binding);
R2S3pw{1,3} = 0;
R2S3pw{2,3} = 0;
R2S3pw{7,3} = 2;
R2S3pw{7,4} = 1;


%%%%cAMP-PKA pathway
Name = [string('Gpr1'),string('Gpa2'),string('Ras1'),string('Ras2'),...
    string('Cdc25'),string('Sdc25'),string('Ira1'),string('Ira2'),...
    string('Cyr1'),string('Tpk1'),string('Tpk2'),string('Tpk3'),string('Bcy1'),...
    string('Rim15'), string('Yak1'),string('Msn2'),string('Msn4'),string('Gis1')]';
presence = ones(length(Name),1);
localization = ones(length(Name),1);
phosphorylation = zeros(length(Name),1);
GXP = zeros(length(Name),1);
DNA_binding  = zeros(length(Name),1);
PKApw= table(Name , presence, localization, phosphorylation, GXP, DNA_binding);
PKApw{1,3} = 0;
PKApw{2,3} = 0;
PKApw{16,3} = 2;
PKApw{17,3} = 2;
PKApw{18,3} = 2;
PKApw{14,3} = 2;
PKApw{14,4} = 0;
PKApw{2,5} = 1;
PKApw{3,5} = 1;
PKApw{4,5} = 1;
placeholders(1) = 0; %placeholder to make model work (0) PKA inactive, not bound to cAMP (1)PKA is bound to cAMP and active

%%%%TOR pathway
Name = [string('Tor1'),string('Tor2'),string('Kog1'),string('Lst8'),string('Tco89'),...
    string('Tap42'),string('PHPs'),string('Sit4'),string('Tip41'),string('Ure2'),...
    string('Gln3'),string('Sch9'), string('Par32'), string('Npr2'), string('Gtr1')]';
presence = ones(length(Name),1);
localization = ones(length(Name),1);
phosphorylation = zeros(length(Name),1);
GXP = zeros(length(Name),1);
DNA_binding = zeros(length(Name),1);
TORpw = table(Name , presence, localization, phosphorylation, GXP, DNA_binding);
TORpw{10,4} = 1; % phosphorelated Ure2
TORpw{11,4} = 1; % phosphorelated Gln3
TORpw{13,4} = 1; % phosphorelated Par32
placeholders(2) = 0; % TORC1

end