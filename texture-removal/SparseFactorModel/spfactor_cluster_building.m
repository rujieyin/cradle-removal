disp('++++++++++++++++++++++++++++++++++++++++++++++++++');
disp(' ');

NumRand = 45;
chunk_size = 1;
work_path = '/home/grad/rachel/matlab/cradleremoval/texture_removal';
search_path = '/home/grad/rachel/matlab/SparseFactorModel';
result_path = '/gtmp/rachel/SparseFactorModel/19a_vertical/';%['/home/grad/rachel/Documents/MATLAB/Preparation2/'];%['/gtmp/rachel/SparseFactorModel/mcmc_spfact/'];
if (~exist(result_path, 'dir'))
    mkdir(result_path);
end
scripts_path = '/gtmp/rachel/SparseFactorModel/script/';%'/home/grad/rachel/Documents/MATLAB/Preparation2/script/';%
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
            fprintf(fid, 'exit; "\n');
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
        fprintf(fid, '#!/bin/bash\n');
        fprintf(fid, '#$ -S /bin/bash\n');
        script_text = ['matlab -nodesktop -nodisplay -nojvm -nosplash -r '...
            '" cd ' work_path ' ;' ...
            ' addpath ' search_path ' ;'];
        fprintf(fid, '%s ',script_text);        
    end
    
    %prepare the script
    
    %using comparison with tps deformation
    script_text = ['BoundaryImg_cluster ' ...
        result_path ' ' ...
        num2str(job_id) '; ' ];
    fprintf(fid, '%s ',script_text);
    
    count=count+1;
    
end


%close the last script file
fprintf(fid, 'exit; "\n');
fclose(fid);

%%qsub
jobname = ['rachel_' num2str(job_id)];
serr = [eo_path 'e_job_' num2str(job_id)];
sout = [eo_path 'o_job_' num2str(job_id)];
tosub = ['!qsub -q all.q -N ' jobname ' -o ' sout ' -e ' serr ' ' script_one_dist_name ];
eval(tosub);

