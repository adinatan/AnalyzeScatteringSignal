function result = sumstructs(fcn, varargin)
% Sum numeric only content of matlab structs with identical internal 
% structure of fields and subfields
%
% example:
% %if the first struct is:
%  f(1).A.a=rand(3);   f(2).A.a=ones(3);
%  f(1).B=1;           f(2).B=2;
% 
%  %and the second struct is:
%  g(1).A.a=rand(3);   g(2).A.a=10*ones(3);
%  g(1).B=10;          g(2).B=20;
% 
%  sum_fg = sumstructs(@(x)sum(x, 3), g, f);
%
% sum_fg(1).B
% ans =
%    11
% sum_fg(2).B
% ans =
%    22
%
%   Ver 1  (2021-12-22)
%   Adi Natan (natan@stanford.edu)

  fcns0 = {@(x)cat(3, x{:}), @(x)x};
  switcher0 = @(tf, s)fcns0{tf+1}(s);
  fcns = {@(s)fcn(s), @(s)sumstructs(fcn, s{:})};
  switcher = @(tf, s)fcns{tf+1}(s);
  c = cellfun(@(x){struct2cell(x)}, varargin);
  s0 = cat(3, c{:});
  s1 = reshape(s0, [], numel(varargin));
  s2 = cellfun(@(x){switcher0(isstruct(x{1}), x)}, num2cell(s1, 2));
  s3 = reshape(s2, size(c{1}));
  s4 = cellfun(@(c){switcher(iscell(c), c)}, s3);
  fnames = fieldnames(varargin{1});
  result = cell2struct(s4, fnames, 1);
end


