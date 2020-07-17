function [smpGirdleOfs] = getGirdleSampleOffset(srng)
%GETGIRDLESAMPLEOFFSET
% Positions in sample girdle are aligned with clock-wise order.
%
pos_size = srng * 8;
edge_len = srng * 2 + 1;
smp_d = srng * 2;
smpGirdleOfs = zeros(pos_size, 4);

% init sampling offsets and grad factor
edge_idx = 1;
for edge_idx = 0:3
    lidx = edge_idx*smp_d + 1;
    ridx = (edge_idx + 1)*smp_d;
    
    if(edge_idx==0)
        smpGirdleOfs(lidx:ridx, 1) = -srng;
        smpGirdleOfs(lidx:ridx, 2) = -srng:1:srng -1;
    elseif(edge_idx==1)
        smpGirdleOfs(lidx:ridx, 1) = -srng:1:srng - 1;
        smpGirdleOfs(lidx:ridx, 2) = srng;
    elseif(edge_idx==2)
        smpGirdleOfs(lidx:ridx, 1) = srng;
        smpGirdleOfs(lidx:ridx, 2) = srng:-1:-srng + 1;
    elseif(edge_idx==3)
        smpGirdleOfs(lidx:ridx, 1) = srng:-1:-srng + 1;
        smpGirdleOfs(lidx:ridx, 2) = -srng;
    end
    
    for i = lidx:ridx
        if(smpGirdleOfs(i, 1)==srng)
            smpGirdleOfs(i, 3) = -1;
        elseif(smpGirdleOfs(i, 1)==-srng)
            smpGirdleOfs(i, 3) = 1;
        end
        
        if(smpGirdleOfs(i, 2)==srng)
            smpGirdleOfs(i, 4) = -1;
        elseif(smpGirdleOfs(i, 2)==-srng)
            smpGirdleOfs(i, 4) = 1;
        end
    end
end

end

