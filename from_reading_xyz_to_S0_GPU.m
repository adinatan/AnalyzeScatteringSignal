clear all
close all

% This part was to read a spesific xyz file... might not be good for you
thepath='';%.\xyz\';
disp('read xyz from file:'); tic;

fnames_in_folder=dir('*.xyz');
for n1=1:numel(fnames_in_folder)
    fn=[thepath  fnames_in_folder(n1).name];
    FileName{n1}= fnames_in_folder(n1).name(1:end-4);
    %disp(['reading file name: ' fn])
    clear c
    fid = fopen(fn);
    n=0;
    while true
        tline = fgetl(fid);
        if ~ischar(tline);  break; end   %end of file
        n=n+1;
        c{n}=strsplit(tline);
        if n<2
            names_length=str2num(c{1}{2});
            Name=cell(1,names_length);
            Namen=zeros(1,names_length,"uint16");
            xyz=zeros(names_length,3,"single");

        end
        if n>2 && n<str2num(c{1}{2})+3
            if numel(c{n})<5
                Name(n-2)=(c{n}(1));
                Namen(n-2)=prod(uint16(c{n}{1}));
                xyz(n-2,:) = cell2mat(cellfun(@str2num,c{n}(2:4),'un',0));

                solute_vec(n-2)=0;
            else
                Name(n-2)=(c{n}(2));
                Namen(n-2)=prod(uint16(c{n}{2}));

                solute_vec(n-2)=1;
                xyz(n-2,:) = single(cell2mat(cellfun(@str2num,c{n}(3:5),'un',0)));
            end
        end
        if n>str2num(c{1}{2})+3 % header size is 2 rows
            break
        end
    end
    % counter ions are deleted
    Name(60:61)=[];
    Namen(60:61)=[];
    xyz(60:61,:)=[];
    fclose('all');
end
disp(['reading done! ' '(' num2str(toc) ' sec)']);


% the assumption from here forward is that we have an array xyz: [# of atoms x 3 dimensions]
xyz=gpuArray(xyz);
% and a vector (Name) of atom names, and\or a vector of numbers (Namen)
% representating an atom names.

%% pair info needs to be done only once: (~ about 6 seconds section)
%disp('set pairs info from atom names:'); tic;
dq=0.01;
q_max=4;
q=eps:dq:q_max;
q=gpuArray(q);

uName =unique(Name);
for n=1:numel(uName);uNamen(n)=prod(single(uName{n}));end

allpairs= gpuArray(nchoosek(uint16(1:numel(Name)),2));

lsa=sum(solute_vec); % last solute atom
% id for solute, solvent and solute_solvent cross terms pairs:
solute_id         =  allpairs(:,1)<=lsa & allpairs(:,2)<=lsa;
solute_solvent_id = (allpairs(:,1)<=lsa & allpairs(:,2)>lsa) | (allpairs(:,1)>lsa & allpairs(:,2)<=lsa) ;
solvent_id        =  allpairs(:,1)>lsa & allpairs(:,2)>lsa;

group_id=[solute_id,solute_solvent_id,solvent_id];
pl= Namen(allpairs);
pl= gpuArray(sort(pl,2));% sort is to make one version of AB and BA pairs
[u_pairs, u_pairs_id, tot_pairs_id] = unique(pl, 'rows');

% form factors function based on the CM coef
f0_CM=@(q,c) sum(c(1:5).*exp(-c(7:11).*(q /(4*pi)).^2))+c(6);

% f0 are the (unique) atomic form factors
for ii=1:numel(uName)
    [cmcoeff, Z]=CMcoef(uName{ii});
    f0(:,ii)=(f0_CM(reshape(q,1,[]),cmcoeff));
end
f0=gpuArray(f0);

% ufafb are the (unique) pair product fa*fb
for ii=1:numel(u_pairs_id)
    ufafb(:,ii)=(f0(:,u_pairs(ii,1)==uNamen).*f0(:,u_pairs(ii,2)==uNamen));
end
ufafb=gpuArray(ufafb); 

%% analysis of pairs to get S0 (need to do every time iter.)
xyz_ab=(xyz(allpairs(:,2),:)-xyz(allpairs(:,1),:));
Rab     = vecnorm(xyz_ab');

% get the memory footprint per pair in a group
for nd=1:3
    for n=1:size(u_pairs,1)
        id=group_id(:,nd)&(tot_pairs_id==n);
        mem_print(n,nd)  =   sum(id)*numel(q) ; % max memory for my gpu is 1.666e+10 bytes
    end
end

% memory chunck to process at at time (in log10)
memory_log_load=8; % this was found staying below "Out of memory" error and optimizing speed

clear pairid groupid
%spair is an array [q x (# of unique pair types) x (# groups: solute, cage, solvent)] capturing the contribution
% to scattering of all pairs per type per group.
spair=zeros(numel(q),size(u_pairs,1),3,"double");  % initialize memory
[pairid, groupid]=find((log10(mem_print).*~isinf(log10(mem_print)))<memory_log_load);

for n=1:numel(pairid)
    id=group_id(:,groupid(n))&(tot_pairs_id==pairid(n));
    ani_order=0;
    qR=q(:)*Rab(id);
    spair(:,pairid(n),groupid(n))=sum(2*(-1)^(ani_order/2)*ufafb(:,pairid(n)).*spherical_bessel_gpu(ani_order,qR ),2);
end

clear pairid groupid
[pairid, groupid]=find(log10(mem_print).*~isinf(log10(mem_print))>=memory_log_load);

for n=1:numel(pairid)
    id=group_id(:,groupid(n))&(tot_pairs_id==pairid(n));
    idx=find(id);
    nn=ceil(sum(id)*numel(q)/(10^memory_log_load));
    idxN=numel(idx);

    splitn = floor(idxN/nn);
    split = 1:splitn:idxN;

    for ni=1:nn-2

        rr=((split(ni):split(ni)+splitn-1));

        ani_order=0;
        qR=q(:)*Rab(idx(rr));
        sum0=gather(sum( 2*(-1)^(ani_order/2)*ufafb(:,pairid(n)).*spherical_bessel_gpu(ani_order,qR ),2 ));
        spair(:,pairid(n),groupid(n))=spair(:,pairid(n),groupid(n))+sum0;

    end
    rr=split(nn-1):idxN;
    qR=q(:)*Rab(idx(rr));

    ani_order=0;
    sum0=gather(sum( 2*(-1)^(ani_order/2)*ufafb(:,pairid(n)).*spherical_bessel_gpu(ani_order,qR ),2 ));
    spair(:,pairid(n),groupid(n))=spair(:,pairid(n),groupid(n))+sum0;


end

%%
% add the atomic part to the scattering (unimportant for the deltaS0 but
% for completness we will do that to get the total isotrpoic scattering)

atom_table=tabulate(Name);
F0=q(:).*0;
for n=1:numel(uName)
    F0=F0+(f0(:,strcmp(atom_table(n,1),uName))).^2.*atom_table{n,2};
end

clear S0
S0=F0+sum(sum(spair,3),2);


disp(['done ! ' '(' num2str(toc) ' sec)']);
semilogy(q,S0);xlabel('Q');ylabel('I [e^2]')
