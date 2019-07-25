function [out1, out2] = cplx_unique(in1,in2,tol)

m = length(in1);
rem_idx = false(m);
for i=1:m-1
    if ~rem_idx(i)
        for j=i+1:m
            if ~rem_idx(j)
                diff = in1(i)-in1(j);
                rem_idx(j) = (abs(real(diff))<tol) & (abs(imag(diff))<tol);
            end
        end
    end
end

out1 = in1;
out2 = in2;

% remove duplicates
out1(rem_idx) = [];
out2(:,rem_idx) = [];

% Sort for real parts
[foo, sortidx] = sort(real(out1),'descend');
out1 = out1(sortidx);
out2 = out2(:,sortidx);

