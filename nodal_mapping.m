% func: mapping layered properties on to grid nodes
% input args:
%               val: vector of layered properties from model
%               nm: vector of numbers of subdivision in each layer
% outputs:
%               nv: vector of nodal values 

function [nv] = nodal_mapping(val,nm)
    if length(val) ~= length(nm)
        error('incorrect nodal mapping input arguments...\n');
    else
        N = sum(nm+1);
        nv = zeros(N,1);
        top = 0;
        for i = 1:length(nm)
            bot = top + 1;
            top = bot + nm(i);
            nv(bot:top) = val(i);
        end
    end   
end