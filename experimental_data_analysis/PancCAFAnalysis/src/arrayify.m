function out = arrayify(T,name,put_field_dims_first)

% take in a structure array and a fieldname of that array and output the
% array of all all elements of the array with their field values arranged
% "nicely". That means, first dimensions will vary over the dimensions of T.
% Then the dimensions of the field will be varied over next.
%
% a variable argument can make the indexes swap order: first vary over
% dimensions of field, then dimensions of T.

arguments
    T
    name
    put_field_dims_first = false
end

if ~isfield(T,name)
    for i = 1:numel(T)
        assert(any(strcmp(fieldnames(T(i)),name))) % for graphics arrays, it doesn't like to collect all the fieldnames together even if they all share one, e.g. YData
    end
end

sz1 = size(T);
sz1(sz1==1) = []; % don't worry about these dimensions

if isempty(sz1) % then T is a single struct
    out = T.(name);
    return;
end

if ~ischar(T(1).(name)) % a single char vector, not a cell array of char vectors
    sz2 = size(T(1).(name));
    sz2(sz2==1) = []; % don't worry about these dimensions either

    if isempty(sz2) % then the fields are scalars
        out = reshape([T.(name)],[sz1,1]);
        return;
    end

    if isnumeric(T(1).(name))
        out = zeros([sz1,sz2]);
    elseif isstring(T(1).(name))
        out = strings([sz1,sz2]);
    elseif all(isgraphics(T(1).(name)))
        out = gobjects([sz1,sz2]);
    end
    sz = size(out);

    out = reshape(out,numel(T),[]);
    for i = 1:numel(T)
        out(i,:) = T(i).(name)(:);
    end

    out = reshape(out,sz);

    if put_field_dims_first==true
        out = permute(out,[(length(sz1)+1):ndims(out),1:length(sz1)]);
    end

else % handle char vectors differently
    out = cell([sz1,1]); % tack on a 1 in case sz1 is a scalar; this prevents matlab from making a square array
    for i = 1:numel(T)
        out{i} = T(i).(name);
    end
end