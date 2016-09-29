% func: create solution file name
% input args: 
%     l0,m0: harmonic of the tide
%     l1,m1: harmonic of the eigenstructure
%     mode,l,m,order: type,degree,order,order of pert of the mode
% outputs:
%     fname: output file name

function [fname] = file_sol(mode,l,m,order,l0,m0,l1,m1)

    global Output;
    dir = Output.dir_sol;
    pre = Output.prefix;

    if mode == 1
        type = 'S';
    elseif mode == -1
        type = 'T';
    end
    fname = sprintf('%s%s.TD%d%d_LH%d%d_%s%d%d_P%d.mat',dir,pre,l0,m0,l1,m1,type,l,m,order);
        
end