clear all
close all
disp('++++++++++++++++++++++++++++++++++++++++++++++++++');
disp(' ');

NumRand = 13;
chunk_size = 1;
work_path = '/home/grad/rachel/matlab/cradleremoval';
search_path = '/home/grad/rachel/matlab/dtcwt';
result_path = ['/home/grad/rachel/Documents/MATLAB/MCAdecomp_Ghent_8_' num2str(NumRand) '/'];
if (~exist(result_path, 'dir'))
    mkdir(result_path);
end
eo_path = [result_path 'err_and_out/'];
if (~exist(eo_path, 'dir'))
    mkdir(eo_path);
end
command_text = ['!rm -f ' eo_path '*']; eval(command_text);
disp(command_text);
scripts_path = [result_path,'script/'];
for k = 1:NumRand
            script_one_dist_name = [scripts_path 'script_' num2str(k)];
    
            jobname = ['racheljob_' num2str(k)];
            serr = [eo_path 'e_job_' num2str(k)];
            sout = [eo_path 'o_job_' num2str(k)];
            tosub = ['!qsub -q all.q -N ' jobname ' -o ' sout ' -e ' serr ' ' script_one_dist_name ];
            eval(tosub);      
        
end

