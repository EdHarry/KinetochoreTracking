function updateSVNrepositories
% This function updates the svn repositories found in %MATLABHOME% and %XTENSIONS%
%
% Aaron Ponti, 2006/11-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get environment variables MATLABHOME (mandatory) and XTENSIONS (optional)
MATLABHOME = getenv( 'MATLABHOME' );
if isempty( MATLABHOME )
    MATLABHOME = getenv( 'HOME' );
    if isempty( MATLABHOME )
        disp( 'Error: environment variable ''MATLABHOME'' or ''HOME'' not set.' );
        return
    end
end

% Check for the existance of the subdirectory 'matlab'
matlabSubdir = fullfile( MATLABHOME, 'matlab' );
if exist( matlabSubdir, 'dir' ) ~= 7
    fprintf(1,'Error: subdirectory %s does not exist.',matlabSubdir);
    return
end

% Get list of subdirectories
subdirs = dir( matlabSubdir );
subdirs( [ subdirs.isdir ] == 0 ) = [];

% Now svn update all directories that are svn working copies
counter = 0;
ccounter = 0;

% check for progressText
try
    progressText(0,'Updating subversion directories') % Create text
    verbose = true;
catch
    verbose = false;
end

nSubdirs = numel( subdirs );
report = cell(nSubdirs,2);
reportCommit = report;
for i = 1 : nSubdirs
    
     
    if strcmp( subdirs( i ).name, '.' ) || strcmp( subdirs( i ).name, '..' ) 
        continue;
    end
    
    % do not check directories to which certain users do not have access
    if strcmp(getenv('COMPUTERNAME'),'MADP05-IRIC') && ~strcmp(getenv('USERNAME'),'p0877743') && ...
        (strcmp(subdirs(i).name,'newFunctions') || strcmp(subdirs(i).name,'chromdyn-tsri'))
        continue
    end

    if verbose
        progressText((i-1)/nSubdirs+eps,sprintf('Updating subversion directories (checking %s)', subdirs( i ).name))
    end
    % Run an svn status first. For some reason, you need to give the full
    % path in order to get the full answer.
    eval(['[ status, err ] = system( ''svn status "',fullfile(matlabSubdir,subdirs( i ).name),'"'' );' ] );

    % Skip directories that are not working copies
    if findstr( err, 'is not a working copy' )
        continue;
    end

    % good directory. Store if something needs committing
    if ~isempty(err)
        ccounter = ccounter + 1;
        reportCommit{ ccounter, 1 } = subdirs( i ).name;
        reportCommit{ ccounter, 2 } = err;
    end
    
    if verbose
        progressText((i-1)/nSubdirs,sprintf('Updating subversion directories (updating %s)', subdirs( i ).name))
    end

    % Run an svn update
    eval(['[ status, err ] = system( ''svn update "',subdirs( i ).name,'"'' );' ] );

    % Store the outputs
    counter = counter + 1;
    report{ counter, 1 } = subdirs( i ).name;
    report{ counter, 2 } = err;

end
if verbose
    progressText(1,'Updating subversion directories');
end

report(counter+1:end,:) = [];
reportCommit(ccounter+1:end,:) = [];


% Save to log

% Directory where to store the file
OUTPUTDIR=MATLABHOME;
% if isempty(OUTPUTDIR)
%     OUTPUTDIR=getenv('TMP');
%     if isempty(OUTPUTDIR)
%         OUTPUTDIR=getenv('TEMP');
%         if isempty(OUTPUTDIR)
%             OUTPUTDIR=getenv('HOME');
%             if isempty(OUTPUTDIR)
%                 OUTPUTDIR=pwd;
%             end
%         end
%     end
% end

% Prepare report. In case multiple users share the same repository,
% everyone has to write their separate logfile, because in Windows, Matlab
% cannot set write permission for all.
tmpfile=fullfile(OUTPUTDIR,sprintf('svn_log_%s.txt',getenv('USERNAME')));

nolog = false;
fid=fopen(tmpfile,'w');
if fid==-1
    nolog = true;
end

% Write to log
if ~nolog
for c = 1 : counter
    fprintf( fid, 'Update of %s:\n %s\n',report{ c, 1 }, report{ c, 2 } );
end
if ccounter > 0
    fprintf( fid, '\n\n Un-committed or locked files:\n');
    for cc = 1 : ccounter;
        fprintf( fid, '%s:\n %s\n',reportCommit{ cc, 1 }, reportCommit{ cc, 2 } );
    end
end

% Close file
fclose( fid );
end

if nolog == false
    if ccounter == 0
        disp( [ 'Done (<a href="file:///',tmpfile,'">details</a>).' ]);
    else
        disp( [ 'Done. Warning: Your working copy contains un-committed changes (<a href="file:///',tmpfile,'">details</a>).' ]);
    end
else
    fprintf('Unable to write file %s\n',tmpfile)
end

