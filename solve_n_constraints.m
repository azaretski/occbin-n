function [sim_lin,sim_pw,sim_ss,oo_ref,M_ref,exitflag]=solve_n_constraints(modnames,rspace,violate,shocks,...
    irfnames,nperiods,curb_retrench,maxiter,init_vars)
% Generalization of a function 'solve_two_constraints' to handle an
% arbitrary number of constraints. The original function developed by Luca
% Guerrieri and Matteo Iacoviello as part of the Occbin package, later
% optimized by the Dynare team.
%
% Suppose there are n regime-switching constraints. Each constraint must
% have two regimes: reference regime (0) and alternative regime (1). The
% inputs and outputs are as follows.
%
% - modnames : either a 2^n-by-1 or a 2-by-1 cell array containing names of
% modfiles corresponding to different regimes, without extensions. In
% either case, modnames{1} must contain the name of the modfile when all
% constraints are in the reference regimes. If modnames is a 2-by-1 cell
% array, then modnames{2} must contain the name of the modfile when only
% the first constraint is in the alternative regime, and all other
% constraints are in the reference regime. Moreover, in this case, the
% names must be formatted as follows. if n = 3 and modnames{1} = 'name',
% then modnames{2} = 'name_100'.
%
% - rspace   : either a 2^n-by-n or an empty [] array with binary values 0
% and 1. Each row corresponds to a specific joint regime, and each colum
% indicates the regime of the corresponding constraint. For example, if n =
% 3, then 000 means that all constraints are in the reference regime, while
% 010 means that only the second constraint is in the alternative regime.
% If both 'modnames' and 'rspace' have 2^n rows, then they must correspond
% to each other: rspace(i,:) defines the regimes present in the model in
% the modfile with name modnames{i}.
%
% Note: the easiest way is to provide a 2-by-1 cell array 'modnames' and an
% empty matrix rspace = []. If the user has formatted the two modfiles in a
% suggested way as explained in the readme file, the function will do
% everything else automatically.
%
% - violate  : 2-by-n cell array of strings (char vectors) that define the
% conditions when the reference regime and the alternative regime are
% violated. A column i corresponds to constraint i, in which case
% violate{1,i} is the violation of the reference regime and violate{2,i} is
% the violation of the alternative regime. See, the Occbin package for more
% details on the required format of the strings.
%
% - exitflag : 1 if the algorithm converged, 0 if not.
%
% The other input and output arguments are the same as in the native Occbin
% functions.
%
% (c) Aliaksandr Zaretski, 2021

global M_ oo_

% Set-up
ncons=size(violate,2);	% number of regime-switching constraints
if isempty(rspace) || size(rspace,1)~=2^ncons || numel(modnames)~=2^ncons
    if numel(modnames)<2
        error('modnames must have at least two elements: reference and first constraint alternative')
    end
    [modnames,rspace]=regimes_setup(ncons,modnames{1},modnames{2});
end
nreg=size(rspace,1);	% total number of regimes

% Run Dynare at each regime
Ms=cell(nreg,1);
oos=cell(nreg,1);
for i=1:nreg
    eval(['dynare ' modnames{i} ' noclearall nolog'])
    oos{i}=oo_;
    Ms{i}=M_;
    if ~isequal(Ms{1}.endo_names,Ms{i}.endo_names)
        error([modnames{1} ' and ' modnames{i} ' need to have identical endogenous variables!'])
    elseif ~isequal(Ms{1}.exo_names,Ms{i}.exo_names)
        error([modnames{1} ' and ' modnames{i} ' need to have identical exogenous variables!'])
    elseif ~isequal(Ms{1}.param_names,Ms{i}.param_names)
        warning([modnames{1} ' and ' modnames{i} ' have different parameter lists!'])
    end
end
oo_ref=oos{1};
M_ref=Ms{1};
nvars=M_ref.endo_nbr;
endog=M_ref.endo_names;
exog=M_ref.exo_names;

% Create reference ss-values and parameters in workspace
sim_ss=oo_ref.dr.ys;
for i=1:nvars
    eval([endog{i} '_ss=sim_ss(i);']);
end
for i=1:M_ref.param_nbr
    eval([M_ref.param_names{i},'= M_ref.params(i);']);
end

% Compute derivatives
cof=NaN(nvars,3*nvars,nreg);
Jbarmat=NaN(nvars,M_ref.exo_nbr,nreg);
Dbarmat=NaN(nvars,nreg);
for j=1:nreg
    Ms{j}.params=M_ref.params;
    [hm1,h,hl1,Jbarmat(:,:,j),residd]=get_deriv(Ms{j},sim_ss);
    cof(:,:,j)=[hm1,h,hl1];
    Dbarmat(:,j)=residd;
end
Dbarmat(:,1)=0;
if isfield(M_ref,'nfwrd')
    [decrulea,decruleb]=get_pq(oo_ref.dr,M_ref.nstatic,M_ref.nfwrd);
else
    [decrulea,decruleb]=get_pq(oo_ref.dr,oo_ref.dr.nstatic,oo_ref.dr.nfwrd);	% older Dynare releases
end

% Process each constraint to uppend a suffix
cons_diff_str=cell(2,ncons);
for j=1:ncons
    cons_diff_str{1,j}=process_constraint(violate{1,j},'_difference',endog,0);
    cons_diff_str{2,j}=process_constraint(violate{2,j},'_difference',endog,0);	% change 0 to 1 to invert inequality
end

% Regime matrix: 0 -- reference, 1 -- alternative
violvecbool=zeros(nperiods+1,ncons);	% rows are time preiods, columns are constraints

% Prepare for loop
if ~exist('curb_retrench','var') || isempty(curb_retrench)
    curb_retrench=0;
end
if ~exist('maxiter','var') || isempty(maxiter)
    maxiter=20;
end
if ~exist('init_vars','var') || isempty(init_vars)
    init_vars=zeros(nvars,1);
end
init_orig=init_vars;
sim_pw=NaN(nperiods,nvars);
wishlist=endog;

% Loop through time periods
exitflag=1;
for t=1:size(shocks,1)
    if mod(t,100)==0
        disp(['t = ' int2str(t)])
    end
    % Iterate until convergence of regimes
    changes=1;
    iter=0;
    while changes && iter<maxiter
        iter=iter+1;
        % Isolate contiguous periods in the other regime
        regimestart=cell(1,ncons);
        for j=1:ncons
            [~,regimestart{j}]=map_regime(violvecbool(:,j));
        end
        % Get data
        sim_lin=mkdatap_anticipated_n(nperiods,decrulea,decruleb,cof,Jbarmat,Dbarmat,regimestart,...
            violvecbool,endog,exog,irfnames,shocks(t,:),init_vars,rspace);
        if isempty(sim_lin)
            exitflag=0;
            break
        end
        % Check regime changes
        for i=1:numel(wishlist)
            eval([wishlist{i}, '_difference=sim_lin(:,i);']);
        end
        newviolvecbool=NaN(nperiods+1,ncons);
        relaxconstraint=NaN(nperiods+1,ncons);
        for j=1:ncons
            newviolvecbool(:,j)=eval(cons_diff_str{1,j});
            relaxconstraint(:,j)=eval(cons_diff_str{2,j});
        end
        if ~any(newviolvecbool(violvecbool==0)) && ~any(relaxconstraint(violvecbool==1))
            changes=0;
        end
        % Update regime matrix
        if curb_retrench
            retrench=zeros(nperiods+1,ncons);   % Gauss-Sidel idea of slowing down the change in the guess
            for j=1:ncons
                retrench(find(relaxconstraint(:,j) & violvecbool(:,j),1,'last'),j)=1;   % only relax last bindig period
            end
            violvecbool=(violvecbool | newviolvecbool)-retrench;
        else
            violvecbool=(violvecbool | newviolvecbool)-(relaxconstraint & violvecbool);
        end
    end
    if ~exitflag
        break
    end
    % Save data and update regimes guess expecting no additional shocks
    init_vars=sim_lin(1,:);
    sim_pw(t,:)=init_vars;
    init_vars=init_vars';
    violvecbool=[violvecbool(2:end,:);zeros(1,ncons)];
end

% Check convergence
if exitflag
    if changes
        disp('Did not converge -- increase maxiter')
        exitflag=0;
    end
    sim_pw(t+1:end,:)=sim_lin(2:nperiods-t+1,:);
end

% Linear solution
sim_lin=mkdata(nperiods,decrulea,decruleb,endog,exog,wishlist,irfnames,shocks,init_orig);

end


% ----------------------------------------------------
% Grid of regimes and modfiles for alternative regimes
% ----------------------------------------------------

function [modnames,rspace]=regimes_setup(ncons,mod_0,mod_1)
% - ncons     : number of regime-switching constraints;
% - mod_0_str : mod file name---without extension---when all
% constraints are in reference regime;
% - mod_1_str : mod file name---without extension---when first
% constraint is in alternative regime and all other constraints are in
% reference regime. This filename must have the following pattern. If n = 3
% and mod_0_str = 'name', then mod_1_str = 'name_100';
% - modnames  : 2^n-by-1 cell vector with modnames for all regimes
% - rspace    : 2^n-by-n grid matrix of regimes for each cosntraint,
% reference regime is 0, alternative regime is 1; row i is the regime of
% the modfile i from modnames

% Construct index string based on number of constraints
ind_1_str=blanks(ncons);
ind_1_str(:)='0';
ind_1_str(1)='1';

% Check that user-provided mod file names are in expected format
if ~strcmp([mod_0 '_' ind_1_str],mod_1)
    error('Mod file names are not in expected format!')
end

% Check that alternative mod file has appropriate regime selection line
reg_1_str=['@#define regimes=[' insert_commas(ind_1_str) ']'];
f=fileread([mod_1 '.mod']);
if ~contains(f,reg_1_str)
    error('Regime selection in alternative mod file is not present or not formatted correctly!')
end

% Create regime space matrix
g=cell(1,ncons);
[g{:}]=ndgrid([0;1]);
rspace=cell2mat(cellfun(@(x)x(:),g,'UniformOutput',false));
nreg=size(rspace,1);	% number of regimes
if nreg==2
    return
end

% Create additional mod files with appropriate regime declaration line; if
% a file exists, check and correct the regime declaration line if necessary
modnames=cell(nreg,1);
modnames{1}=mod_0;
modnames{2}=mod_1;
for i=3:nreg
    ind_i_str=sprintf('%d',rspace(i,:));
    modnames{i}=[mod_0 '_' ind_i_str];
    fname=[modnames{i} '.mod'];
    if ~exist(fname,'file')
        copyfile([mod_1 '.mod'],fname)
    end
    reg_i_str=['@#define regimes=[' insert_commas(ind_i_str) ']'];
    f=fileread(fname);
    if ~contains(f,reg_i_str)
        f=replace(f,reg_1_str,reg_i_str);   % appropriate regime
        x=fopen(fname,'w');
        fprintf(x,'%s',f);
        fclose(x);
    end
end

end


% -----------------------------------------------------------------
% Function that inserts a commma between each character of a string
% -----------------------------------------------------------------

function sc=insert_commas(s)
% - s  : character vector, e.g., '123'
% - sc : corresponding vector with commas, e.g., '1,2,3'

if ~ischar(s)
    error('Not character input!')
end

ns=numel(s);
if ns<=1
    sc=s;
else
    sc=blanks(2*ns-1);
    p=1;
    for i=1:ns-1
        sc(p:p+1)=[s(i) ','];
        p=p+2;
    end
    sc(end)=s(end);
end

end
