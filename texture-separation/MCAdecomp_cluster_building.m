clear all;
close all;

disp('++++++++++++++++++++++++++++++++++++++++++++++++++');
disp(' ');

NumRand = 8;
chunk_size = 1;
SampNumb = [400, 400];
work_path = '/home/grad/rachel/matlab/cradleremoval/texture_separation';
search_path = '/home/grad/rachel/matlab/dtcwt';
% % have all maximal conformal factors in sample
result_path = ['/home/grad/rachel/Documents/MATLAB/MCAdecomp_Ghent_8/'];
if (~exist(result_path, 'dir'))
    mkdir(result_path);
end

scripts_path = [result_path 'scripts/'];
if (~exist(scripts_path, 'dir'))
    mkdir(scripts_path);
end

eo_path = [result_path 'err_and_out/'];
if (~exist(eo_path, 'dir'))
    mkdir(eo_path);
end

command_text = ['!rm -f ' scripts_path '*']; eval(command_text);
disp(command_text);
command_text = ['!rm -f ' eo_path '*']; eval(command_text);
disp(command_text);

count = 0;
job_id = 0;
for k = 1:NumRand
    if( mod(count,chunk_size) == 0)  %save the old script and create new script
        
        if (job_id > 0 ) %not the first time
            %close the script file (except the last one, see below)
            fprintf(fid, '%s\n','exit; " ');
            fid
            fclose(fid);
            
            %%qsub
            jobname = ['rachel_' num2str(job_id)];
            serr = [eo_path 'e_job_' num2str(job_id)];
            sout = [eo_path 'o_job_' num2str(job_id)];
            tosub = ['!qsub -q all.q -N ' jobname ' -o ' sout ' -e ' serr ' ' script_one_dist_name ];
            eval(tosub);
        end
        
        job_id = job_id + 1;
        script_one_dist_name = [scripts_path 'script_' num2str(job_id)];
        
        %open the next (first?) script file
        fid = fopen(script_one_dist_name,'w');
        fprintf(fid,'%s\n', '#!/bin/bash');
        fprintf(fid, '%s\n','#$ -S /bin/bash');
        script_text = ['matlab -nodesktop -nodisplay -nojvm -nosplash -r '...
            '" cd ' work_path ' ;' ...
            ' addpath ' search_path ' ;'];
        fprintf(fid, '%s ',script_text);
        
    end
    
    %prepare the script
    
    script_text = [' MCAdecomposition_cluster ' ...
        num2str(k) ' ' ...
        result_path ' ' ...
        num2str(job_id) '; ' ];
    fprintf(fid, '%s ',script_text);
    
    count=count+1;
    
end


%close the last script file
fprintf(fid,'%s\n', 'exit; " ');
fclose(fid);

%%qsub
jobname = ['rachel_' num2str(job_id)];
serr = [eo_path 'e_job_' num2str(job_id)];
sout = [eo_path 'o_job_' num2str(job_id)];
tosub = ['!qsub -q all.q -N ' jobname ' -o ' sout ' -e ' serr ' ' script_one_dist_name ];
eval(tosub);




