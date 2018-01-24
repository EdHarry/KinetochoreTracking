function practiceDataArray
% practice with cells, structure, and dataset arrays

% see Statistics Toolbox / Organizing Data for info on dataset arrays

% 2009 Kathryn Applegate

%% create cell array of lab members
labCell = {...
'Gaudenz'       'Danuser'       'gdanuser@scripps.edu'  'M' 'Principal Investigator';...
'Khuloud'       'Jaqaman'       'kjaqaman@scripps.edu'	'F'	'Research Associate';...
'Dinah'         'Loerke'        'dloerke@scripps.edu'	'F'	'Research Associate';...
'Alex'          'Matov'         'amatov@scripps.edu'	'M'	'Graduate Student';...
'Cristen'       'Kern'          'cristen@scripps.edu'	'F'	'Research Assistant';...
'James'         'Lim'           'jilim@scripps.edu'     'M'	'Graduate Student';
'Mohsen'        'Sabouri Ghomi' 'sabouri@scripps.edu'	'M'	'Research Associate';...
'Andrea'        'Bacconi'       'abacconi@scripps.edu'	'M'	'Research Associate';...
'Hunter'        'Elliott'       'helliott@scripps.edu'	'M'	'Graduate Student';...
'Kathryn'       'Applegate'     'applekat@scripps.edu'	'F'	'Graduate Student';...
'Danny'         'Nunez'         'dnunez@ucsd.edu'       'M'	'Graduate Student';...
'Pei-Hsin'      'Hsu'           'peihsin@scripps.edu'	'M'	'Graduate Student';...
'Shann-Ching'	'Chen'          'shanncc@scripps.edu'	'M'	'Research Associate';...
'Sylvain'       'Berlemont'     'sylvain@scripps.edu'	'M'	'Research Associate';...
'Maria'         'Bagonis'       'mbagonis5@gmail.com'	'F'	'Graduate Student';...
'Allen'         'Liu'           'allenliu@scripps.edu'	'M'	'Research Associate';...
'Miriam'		'Berba'         'mirberba@scripps.edu'	'F'	'Administrative Assistant';...
'Zion'          'Maffeo'        'zmaffeo@scripps.edu'	'M'	'Systems Administrator'};

nMembers = size(labCell,1);

% find subset - all the graduate students (Method 1)
temp=zeros(nMembers,1);
for iMember=1:nMembers
    if strcmp(labCell{iMember,5},'Graduate Student')
        temp(iMember)=1;
    end
end
matchIdx=find(temp);
gradStudents=labCell(matchIdx,:);

% find subset - all the postdocs (Method 2)
matchIdx1=cellfun(@(x) strcmp(x,'Research Associate'),labCell(:,5));
matchIdx2=arrayfun(@(x) strcmp(x,'Research Associate'),labCell(:,5));
postDocs=labCell(matchIdx1,:);


%% convert cell array to structure format
labStruct=cell2struct(labCell,{'firstName','lastName','email','sex','jobDescription'},2);

% find subset - all the females (Method 1)
temp=zeros(nMembers,1);
for iMember=1:nMembers
    if strcmp(labStruct(iMember).sex,'F')
        temp(iMember)=1;
    end
end
matchIdx=find(temp);
females=labStruct(matchIdx);

% find subset - all the postdocs (Method 2)
matchIdx=arrayfun(@(x) strcmp(x.jobDescription,'Research Associate'), labStruct);
postDocs=labStruct(matchIdx);


%% sorting

% sort cell array by last names
[temp,idx]=sort([labCell(:,2)]);
labCellOrdered=labCell(idx,:);

% sort structure by last names
lastNameCell={labStruct.lastName}';
[temp,idx]=sort(lastNameCell);
labStructOrdered=labStruct(idx);


%% create dataset array from structure

% re-order the job description labels with categorical arrays
job = labCell(:,5);
job1 = ordinal(job);
jobLabels1 = getlabels(job1)';
job2 = reorderlevels(job1,jobLabels1([3 5 2 4 6 1]));
jobLabels2 = getlabels(job2)';

lastName=nominal(labCell(:,2));
nameLabels = getlabels(lastName);

% make label for each member
memberNum = strcat({'Member'},num2str((1:nMembers)','%d'));

% create the dataset array
labDataset = dataset({labCell(:,1:4),'firstName','lastName','email','sex'},{job2,'jobDescription'},...
    'ObsNames',memberNum);

% set some info and then retrieve it
desc = 'Info about LCCB members';
info = 'http://lccb.scripps.edu/';
dims = {'lab members','misc info'}
labDataset = set(labDataset, 'Description',desc, 'UserData',info, 'DimNames',dims);

get(labDataset)

% sort array easily
alphabetical = sortrows(labDataset,'lastName');        
jobRanking = sortrows(labDataset,'jobDescription');


% easy access for subsets
labDataset.email(1:5)
labDataset({'Member 1', 'Member 3', 'Member 5'},'lastName')

% joining dataset arrays
a = labDataset(1:5,2);
b = labDataset(1:5,5);
c = [b a];

% broadcast info about chromosomes to new dataset
snames = {'F';'M'};
chrom = dataset({snames,'sex'},{['xx';'xy'],'chromosomes'});
chromDataset = join(labDataset(:,1:3:4),chrom);


