function varargout = make_dep(depstr)
% MAKE_DEP
% Construct dependent variables (e.g. x, vx, ekin, etc) from string input 
% provided by smilei

dep1_str = depstr;
dep1_tok = split(dep1_str);
if str2num(dep1_tok{5}) == 1
  dep1 = logspace(str2num(dep1_tok{2}),str2num(dep1_tok{3}),str2num(dep1_tok{4}));
else
  dep1 = linspace(str2num(dep1_tok{2}),str2num(dep1_tok{3}),str2num(dep1_tok{4}));
end

if nargout == 1  
  varargout{1} = dep1;
elseif nargout == 2
  varargout{1} = dep1_tok{1};
  varargout{2} = dep1;
end