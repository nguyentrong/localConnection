function N = numContactHypercolumn(preType,postType,synDist)
    % input : preType,postType,synDist
    % where preType, postType = type of pre/postsynaptic neuron
    % synDist = distance between 2 synapse to form a synaptic connection
    % (default synDist = 1um see Hellwig 1998)
    % output: N
    % Estimation the number of synaptic contacts between a hypercolumn (size 100umx100um)
    % and its surroundings in 2D plane sing knnsearch 
    % example: 
    % Find 2 nearest neighbors in X and the corresponding values to each
    % point in Y using the distance metric 'cityblock'
    % X = randn(100,5);
    % Y = randn(25, 5);
    % [idx, dist] = knnsearch(X,Y,'dist','cityblock');
          
    global dimx dimy dimz
   
    % size of the hypercolumn in um
    dimx = 200;   
    dimy = 200;   
    dimz = 200;   

    % Neuron density (Braitenberg&Schuez 1998)
    nDensity = 3e4; % neurons/mm^3

    % Number of neurons / hypercolumn
    fprintf('# Estimating the number of neurons...');
    nNumber = nDensity * dimx/1e3 * dimy/1e3 * dimz/1e3;
    pyNum = floor(0.85*nNumber); % Number of pyramidal neurons = 85% total neuron number
    inNum = floor(0.15*nNumber/2); % The rest are interneurons
    if strcmp(preType,'P'); preNum = pyNum; else preNum = inNum;  end  
    if strcmp(postType,'P'); postNum = pyNum; else  postNum = inNum; end
    fprintf('completed\n');    
    
    %Support that the somas normally distribute  in the hypercolumn 
    % -> Find the location of each soma
    fprintf('# Estimating locations of the somas...');    
    preSoma = locationSoma(preNum);
    postSoma = locationSoma(postNum);
    fprintf('completed\n');
    
    % Load neurons into the hypercolumn
    fprintf('# Loading neurons into the hypercolumn...');    
    preNeuronData  = loadNeuron(preType,'A',preSoma);
    postNeuronData = loadNeuron(postType,'D',postSoma);
    fprintf('completed\n');
    
    % Now estimate the number of synaptic contact depending on the
    % separation between the pre and postsynaptic hypercolumn. Note that
    % the range of connection wont be greater then 500um, here use this
    % limit 1000um = 5x hypercolumn size 
    N = [];
    fprintf('# Prosessing (total 121):     \t\t');
    counter = 1;
    for sepx = -5:5
       for sepy = -5:5
            x = sepx*dimx; y = sepy*dimy;
            dat = [postNeuronData(:,1)+x postNeuronData(:,2)+y postNeuronData(:,3)];
            [~, dist] = knnsearch(preNeuronData,dat,'dist','euclidean');
            dist = dist(dist<=synDist);
            if size(dist,1) > 0
                N = [N; sepx sepy size(dist,1)]; 
            end
            fprintf('\b\b\b%03d',counter);
            counter = counter + 1;
       end
    end    
end

function SOMA = locationSoma(num)
    global dimx dimy dimz
    n = dimx*dimy*dimz / num;
    d = n^(1/3);
    SOMA = [];
    for idx = 1:ceil(dimx/d)
        for idy = 1:ceil(dimy/d)
            for idz = 1:ceil(dimz/d)
                SOMA = [SOMA; idx*ceil(dimx/d) idy*ceil(dimy/d) idz*ceil(dimz/d)];
            end
        end
    end
end

function NEURON = loadNeuron(type,c,soma_location)
    d = dir([type '/' c]);
    nfiles = size(d,1) - 2;
    n = size(soma_location,1);
    % choose random file
    fid = 1 + floor((nfiles-1).*rand(n,1));
    NEURON = [];
    for i = 1:n
        somx = soma_location(i,1);
        somy = soma_location(i,2);
        somz = soma_location(i,3);
        dat = importdata([type '/' c '/' c '' num2str(fid(i)) '.txt']);
        if ~isempty(dat)
            NEURON = [NEURON; dat(:,1)+somx dat(:,2)+somy dat(:,3)+somz];       
        end
    end
end































