function sim_paths=mkdatap_anticipated_n(nperiods,decrulea,decruleb,cof,Jbarmat,Dbarmat,regimestart,...
    violvecbool,endog,exog,irfnames,scalefactor,init_vars,rspace)
% Generalization of a function 'mkdatap_anticipated_2constraints' to handle
% an arbitrary number of constraints. The original function developed by
% Luca Guerrieri and Matteo Iacoviello as part of the Occbin package, later
% optimized by the Dynare team.
%
% This function is called exclusively by another function
% 'solve_n_constraints' and is provided as a separate file only for the
% sake of readability. For the nature of inputs and outputs, see
% 'solve_n_constraints' and the original Occbin analogs.
%
% (c) Aliaksandr Zaretski, 2021

% Set-up
nreg=size(rspace,1);	% total number of regimes
ncons=size(rspace,2);	% number of regime-switching constraints
nvars=size(endog,1);
if ~exist('scalefactor','var') || isempty(scalefactor)
    scalefactor=ones(size(irfnames,1),1);
end
if ~exist('init_vars','var') || isempty(init_vars)
    init_vars=zeros(nvars,1);
end

% Create matrices for each regime
Cbarmat=NaN(nvars,nvars,nreg);
Bbarmat=NaN(nvars,nvars,nreg);
Abarmat=NaN(nvars,nvars,nreg);
for i=1:nreg
    Cbarmat(:,:,i)=cof(:,1:nvars,i);
    Bbarmat(:,:,i)=cof(:,nvars+1:2*nvars,i);
    Abarmat(:,:,i)=cof(:,2*nvars+1:3*nvars,i);
end

% Position of the last period when any constraint binds
last_regime_start=NaN(1,ncons);
for i=1:ncons
    last_regime_start(i)=regimestart{i}(end);
end
Tmax=max(last_regime_start)-1;

% Get the time-dependent decision rules
if Tmax>0
    indr=NaN(Tmax,1);	% regime indices
    for t=1:Tmax
        indr(t)=find(sum(violvecbool(t,:)==rspace,2)==ncons);
    end
    P=zeros(nvars,nvars,Tmax);
    D=zeros(nvars,Tmax);
    lastwarn('','');	% reset warnings
    invmat=inv(Abarmat(:,:,indr(Tmax))*decrulea+Bbarmat(:,:,indr(Tmax)));
    if ~isempty(lastwarn)
        sim_paths=[];
        return
    end
    P(:,:,Tmax)=-invmat*Cbarmat(:,:,indr(Tmax));
    D(:,Tmax)=-invmat*Dbarmat(:,indr(Tmax));
    if Tmax>1
        for t=Tmax-1:-1:1
            invmat=inv(Bbarmat(:,:,indr(t))+Abarmat(:,:,indr(t))*P(:,:,t+1));
            if ~isempty(lastwarn)
                sim_paths=[];
                return
            end
            P(:,:,t)=-invmat*Cbarmat(:,:,indr(t));
            D(:,t)=-invmat*(Abarmat(:,:,indr(t))*D(:,t+1)+Dbarmat(:,indr(t)));
        end
    end
    E=-invmat*Jbarmat(:,:,indr(1));	% using last computed invmat
end

% Generate data
history=NaN(nvars,nperiods+1);    % state vectors stored columnwise
history(:,1)=init_vars;
errvec=zeros(numel(exog),1);
for i=1:size(irfnames,1)
    shockpos=strcmp(deblank(irfnames(i,:)),exog);
    if any(shockpos)
        errvec(shockpos)=scalefactor(i);
    else
        error(['Shock ',irfnames(i,:),' is not in the model!']);
    end
end

% Deal with shocks
if Tmax>0
    history(:,2)=P(:,:,1)*history(:,1)+D(:,1)+E*errvec;
else
    history(:,2)=decrulea*history(:,1)+decruleb*errvec;
end
if Tmax>1
    for t=2:Tmax
        history(:,t+1)=P(:,:,t)*history(:,t)+D(:,t);
    end
end
for t=max(2,Tmax+1):nperiods+1
    history(:,t+1)=decrulea*history(:,t);
end
history=history';
sim_paths=history(2:end,:);

end
