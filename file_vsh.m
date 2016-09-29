% func: create output file name for vsh_expan
% input args: 
%     l0,m0: harmonic of the "parent" mode
%     l1,m1: harmonic of the eigenstructure
%     mode:  spheroidal (1) or toroidal (-1) of the "parent" mode
%     dir: output file directory
% outputs:
%     fname: output file name

function [fname] = file_vsh(l0,m0,l1,m1,mode,dir)

    if ~strcmp(dir,'./vsh_expan/')
        dir = './vsh_expan/';
    end   
    if mode == 1
        type = 'S';
    elseif mode == -1
        type = 'T';
    end
    fname = sprintf('%s%s%d%d_LH%d%d.mat',dir,type,l0,m0,l1,m1);
        
end