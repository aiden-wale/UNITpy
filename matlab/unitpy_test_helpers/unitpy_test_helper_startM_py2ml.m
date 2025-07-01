function M = unitpy_test_helper_startM_py2ml(m)


if isfield(m, 'in')
    nin  = length(m.in);

    M = m;
    M = rmfield(M, 'in');
    % M = rmfield(M, 'out');

    for k=1:nin
        M.in(k) = m.in{k};
    end
end

% if isfield(m, 'out')
%     nout = length(m.out);

%     M = m;
%     M = rmfield(M, 'out');
%     % M = rmfield(M, 'out');

%     for k=1:nout
%         M.out(k) = m.out{k};
%     end
% end

if isfield(m, 'out')
    fn = fieldnames(m.out);
    for i=1:length(fn)
        M.out = setfield(M.out, fn{i}, getfield(m.out, fn{i}));
    end
end

M = orderfields(M);
