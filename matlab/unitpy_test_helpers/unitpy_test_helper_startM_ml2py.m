function M = unitpy_test_helper_startM_ml2py(m)


if isfield(m, 'in')
    nin  = length(m.in);

    M = m;
    M = rmfield(M, 'in');
    % M = rmfield(M, 'out');

    M.in = cell(1,nin);
    for k=1:nin
        fn = fieldnames(m.in(k));
        for i=1:length(fn)
            M.in{k} = setfield(M.in{k}, fn{i}, getfield(m.in(k), fn{i}));
        end
    end
end

% if isfield(m, 'out')
%     nout = length(m.out);

%     M = m;
%     M = rmfield(M, 'out');
%     % M = rmfield(M, 'out');

%     M.out = cell(1,nout);
%     for k=1:nout
%         fn = fieldnames(m.out(k));
%         for i=1:length(fn)
%             M.out{k} = setfield(M.out{k}, fn{i}, getfield(m.out(k), fn{i}));
%         end
%     end
% end

if isfield(m, 'out')
    fn = fieldnames(m.out);
    for i=1:length(fn)
        M.out = setfield(M.out, fn{i}, getfield(m.out, fn{i}));
    end
end

M = orderfields(M);
