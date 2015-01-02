clc; clear;
slices = 25;
symbols = ['a':'z' 'A':'Z' '0':'9'];
rng('shuffle');

fname = symbols(randi(numel(symbols),[1 8]));

% Also wider parameter set
minn = log10([6,1E-20,1E-5,1E-20,1E-5,... % U2,xFwd1,xRev1,xFwd3,xRev3
    1E-3,1E-3,1E-4,1E-4,1E-1,10,1E-6, 1]); % AXLint1,AXLint2,kRec,kDeg,fElse,AXL,Gas

maxx = log10([1E5,1E5,1E5,1E5,1E5,... % U2,xFwd1,xRev1,xFwd3,xRev3
    1,1,0.1,1,1,1E5,1,10]); % 'AXLint1,AXLint2,kRec,kDeg,fElse,AXL,Gas

minn = [minn minn(1:5)];
maxx = [maxx maxx(1:5)];

Dopts = psoptimset('TimeLimit',60*60,'Display','final',...
    'UseParallel',true,'CompletePoll', 'on');

A = zeros(4,length(minn));
A(1,2:4) = [1, 0, -1];
A(2,3:5) = [-1, 0, 1];
A(1,(2:4) + 13) = [1, 0, -1];
A(2,(3:5) + 13) = [-1, 0, 1];
b = [0; 0; 0; 0];

parpool(4);

for xxxx = 1:100
    for ii = 0:(slices*length(maxx) - 1)
        IDX = mod(ii,length(maxx))+1;
        vv = linspace(minn(IDX),maxx(IDX),slices);

        vIDX = floor(ii/length(maxx)) + 1;
        
        disp([IDX vIDX]);

        minn2 = minn; minn2(IDX) = vv(vIDX);
        maxx2 = maxx; maxx2(IDX) = vv(vIDX); 

        params = minn2 + (rand(size(minn2)) .* (maxx2 - minn2));
        params(4) = params(2);
        params(5) = params(3);
        params(4+13) = params(2+13);
        params(5+13) = params(3+13);
        
        try
            [paramOpt(ii+1,:),fitIDXglobal(ii+1)] = ...
                patternsearch(@cLibA549both,params,A,b,[],[],minn2,maxx2,[],Dopts);
        catch
           paramOpt(ii+1,:) = params;
           fitIDXglobal(ii+1) = 1E6;
        end
    end
    
    fitStruct{xxxx}.paramOpt = paramOpt;
    fitStruct{xxxx}.fitIDXglobal = fitIDXglobal;

    save(['both' fname]);
end