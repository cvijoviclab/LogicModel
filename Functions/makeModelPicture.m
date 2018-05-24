function makeModelPicture(path, dataName)

% creates a plot of the cell with proteins that are defined in the data located
% at 'path' with name 'dataName'

% matrix which consists of the columns: 
% 'name': name of the protein
% 'presence': (1) presencent, (0) not presencent
% 'localization': localizationation of the protein (2) in the nucleus, (1) in the cytosol, 
%           (0) membranebound
% 'phos': (1) phosphorenlated, (0) not phosphorelated
% 'GXP': (2) GTP, (1) GDP, (0) none of them
% 'DNA_binding': (1) yes, (0) no

    pathToData = [path, dataName];

    % load data
    dataTable = readtable(pathToData,'Delimiter','\t');
    data = table2struct(dataTable, 'ToScalar', true);
    tag = numel(fieldnames(data)) > 3; % defines if tags are actually included in data

    % plot cell
    plotEmptyCell();

    % define coordinates for certain positions in the cell
    nProteinsMembranebound = sum(data.presence == 1 & data.localization == 0);
    nProteinsInCytosol = sum(data.presence == 1 & data.localization == 1);
    if tag == 1
        nProteinsNucleus = sum(((data.presence == 1 & data.localization == 2) + (data.DNA_binding == 0)) == 2);
        nProteinsNucleusDNA_binding = sum(((data.localization == 2) + (data.DNA_binding == 1)) == 2);
    else
        nProteinsNucleus = sum(data.presence == 1 & data.localization == 2);
        nProteinsNucleusDNA_binding = 0;
    end

    [cInNucleus, cInNucleusDNA_binding, cInCytosol, cMembrane] = ...
        getProteinPositions(nProteinsNucleus, nProteinsNucleusDNA_binding, nProteinsInCytosol, ...
        nProteinsMembranebound);

    % include proteins
    nProteins = length(data.Name);
    count = [1 1 1 1];

    for itP = 1:nProteins
        proteinName = data.Name{itP};

        if data.presence(itP) >= 1 % only go on if protein is present

            % check localizationation 
            % either in (1) nucleus, (2) in nucleus and DNA bound, (3) in cytosol or
            % (4) membrane bound
            if data.localization(itP) == 2 && data.DNA_binding(itP) == 0
                position = cInNucleus(:,count(1));
                count(1) = count(1) + 1;
            elseif data.localization(itP) == 2 && tag == 1 && data.DNA_binding(itP) == 1
                position = cInNucleusDNA_binding(:,count(2));
                count(2) = count(2) + 1;
            elseif data.localization(itP) == 1
                position = cInCytosol(:,count(3));
                count(3) = count(3) + 1;
            else 
                position = cMembrane(:,count(4));
                count(4) = count(4) + 1;
            end

            % check if protein is modified (phosphorelated, GTP, GDP and DNA bound)
            if tag == 1 
                [tagName, DNA_binding] = getTagName(data, itP);
            else 
                tagName = '';
                DNA_binding = false;
            end

            % place protein in plot
            includeProtein(position, proteinName, tagName, DNA_binding);

        end  

    end

    % save picture
    name = dataName(1:end-4); 
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto', 'PaperSize',[pos(3),pos(4)]);
%     title(['\textbf{',name,'}'], 'FontSize', 32, 'Interpreter' , 'latex')
%     set(get(gca,'title'),'Position',[-1.2 3.45])
    print(gcf,[path,'Pictures/', name], '-depsc', '-r0')
    close(gcf)
    
end