%%% File name:              CHDsolver_TravellingWaves_LiebLinigerGases.m (a Matlab function)
%%% Description:            ... 
%%% Principal developer:    Seyed Nariman Saadatmand 
%%% Contact:                n.saadatmand@griffith.edu.au
%%% Created in:             2/Aug/2018


function CHDsolver(gamma_bg,beta,sigma_bar,N_bg,dt,N_t)

    %%% Function inputs description:
    % [...]		...
    
    %%% Main output file details:
    % [...]		...
    
    
    %%% Setting some parameters for global usage:
    DIR='/home/nariman/Dropbox/AcademiaJobs-eDesktop/MyPapers/1dBoseGases-Qshockwaves/MatlabCollection';
    %DIR='.';  
     
    %N_save = 100;                      % at how many time step, the program should output a data file.
    L=1.0;                              % the size of box (we prefer to set L as the unit of length).
    ParticleDens_ColorMap = 'no';       % whether or not to print a colormap data file.        
    %setenv('EDITOR','vim');            % setting the text editor for external runs.
    N_x = 500;                          % number of spatial discretization points (it is desirable to keep this as par with parallel tDMRG calculations, if existing).
    dx = L/(N_x-1);                     % spatial infinitesimal size used in discretization.  
    x_line = (-L/2:dx:L/2);             % a handy vector that represents spatial discretization.
    tol_iter = 1e-14;                   % tolerance for iterative fix-point solver method(s).
    max_iter = 1000;                    % maximum iteration number to stop fix-point solvers upon not converging.
    rho = zeros(1,N_x);                 % initializing temparory density array.
    density_array = zeros(N_x,N_t+1);   % initializing full density array.
    t_array = zeros(N_t+1,1);           % initializing time-array in our units.
    x_array = zeros(N_x,1);             % initializing position-array in our units.


    if strcmp(ParticleDens_ColorMap,'yes')
        filename_ParticleDens_ColorMap = strcat('ParticleDens_ColorMap-LiebLiniger-GammaBG_',num2str(gamma_bg),'-dt_',num2str(dt),'-CHDsolver_TravellingWaves.out');          
        if exist(fullfile(DIR,filename_ParticleDens_ColorMap),'file') 
          fprintf('NOTE: file %s already exist; if any current data are available, will be used for calculations ...\n', fullfile(DIR,filename_ParticleDens_ColorMap));
          FileID_ParticleDens_ColorMap = fopen( fullfile(DIR,filename_ParticleDens_ColorMap) , 'at');
          if FileID_ParticleDens_ColorMap==-1
            error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,filename_ParticleDens_ColorMap));
          end
        else  
          edit(fullfile(DIR,filename_ParticleDens_ColorMap));
          FileID_ParticleDens_ColorMap = fopen( fullfile(DIR,filename_ParticleDens_ColorMap) , 'at');
          if FileID_ParticleDens_ColorMap==-1
            error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,filename_ParticleDens_ColorMap));
          end
          fprintf(FileID_ParticleDens_ColorMap,'#x/L (our unit)\t#time (Damski''s unit)\t#rho*L (our unit)\n'); % printing the file header
        end 
    end
    
    
    %%% Setting up all required initial conditions for tt=0:                
    %sigma_bar = 0.059;
    %beta = 0.5*(r-1)/(1-r*exp(-1.0/(8*sigma^2)));
    %u0 = (1-rho0*N_x) / sum( exp(- x_line.^2 / (2*kappa^2) ) );             
    %g_bar = 7.5e3;
    g_bar = gamma_bg*N_bg^2*(L+sqrt(2*pi)*beta*(sigma_bar*L)*erf(1/(sqrt(8.0)*sigma_bar*L)));
    %N_bgD = 1 / ( L + beta*sqrt(2*pi)*(sigma_bar*L) );
    N_bgD = 1/(L+sqrt(2*pi)*beta*(sigma_bar*L)*erf(1/(sqrt(8.0)*sigma_bar*L)));
    N = N_bg*(L+sqrt(2*pi)*beta*(sigma_bar*L)*erf(1/(sqrt(8.0)*sigma_bar*L)));
    %fprintf('NOTE: calculating for g_bar=%.4g, N_bg=%.4g, beta=%.4g, sigma_bar=%.4g, and gamma_bg=%.4g ...\n', g_bar, N_bg, beta, sigma_bar, gamma_bg); %DEBUG
    
%     tDMRG_SourceDir='~/Dropbox/AcademiaJobs-eDesktop/MyPapers/1dBoseGases-Qshockwaves/DataCollection_mirrored/ParticleDensity/BoseHubbardU1'; % setting the source directory.
%     FilenameCurrent ='Hdimple_quench-BoseGases-L200q30-alpha1_3.5-m100-BoseHubbardU1-U_0.1-V0_0.05.psi.t0';
%     fileID_tDMRG = fopen(fullfile(tDMRG_SourceDir,FilenameCurrent),'r');
%     formatSpec = '%i %f %f %f %f';
%     sizeIN = [5 inf];
%     INarray = fscanf(fileID_tDMRG,formatSpec,sizeIN);
%     fclose(fileID_tDMRG);
    for xx = 1:N_x
       rho(xx) = N_bgD * (1 + beta*exp(-x_line(xx)^2/(2*(sigma_bar*L)^2)));
       %V0 = (V0_in-V0_out)*exp(-x_line(xx)^2/BoxHalfWidth^2) + V0_out;  
    end
%     N=30;
%     rho(:) = INarray(4,:)/N;
      
    %%% preparing and print out the initial results:
%    filename0 = strcat('rho_x-LiebLiniger-GammaBG_',num2str(gamma_bg),'-dt_',num2str(dt),'-CHDsolver_TravellingWaves.t','0');
%    if exist(fullfile(DIR,filename0),'file')
%     fprintf('NOTE: file %s already exist; new data will be attached to its end ...\n', fullfile(DIR,filename0));
%     FileID0 = fopen( fullfile(DIR,filename0) , 'wt');
%     if FileID0==-1
%       error('ERROR: cannot open the following file for writing, ''%s''.', fullfile(DIR,filename0));
%     end
%    else
%     edit(fullfile(DIR,filename0));
%     FileID0 = fopen( fullfile(DIR,filename0) , 'at');
%     if FileID0==-1
%       error('ERROR: cannot open the following file for writing, ''%s''.', fullfile(DIR,filename0));
%     end
%     fprintf(FileID0,'#x/L(our unit)\t#rho*L(our unit)\n');   % printing the file header 
%    end 
    
    %disp(size(x_line)); % DEBUGGING
    %disp(size(rho)); % DEBUGGING
    %disp(size(rho_particles)); % DEBUGGING
    
    fprintf('NOTE: printing out the initial rho ...\n'); % all will be saved as a .mat file.
    for xx = 1:N_x
       %fprintf('%.16f\t%.16f\n', x_line(xx), rho_particles(xx));
       %fprintf(FileID0, '%.16f\t%.16f\n', x_line(xx)/L, rho(xx)*L*N_bg*(1+sqrt(2*pi)*beta*sigma_bar*erf(1/(sqrt(8.0)*sigma_bar))));
       x_array(xx) = x_line(xx)/L;
       density_array(xx,1) = rho(xx)*N; 
       %save('full_rho.mat','full_rho');
    end 
%    fclose(FileID0);
    
    if strcmp(ParticleDens_ColorMap,'yes')
      for xx = 1:N_x
         %fprintf('%.16f\t%.16f\t%.16f\n', x_line(xx), 0, rho_particles(xx));
         fprintf(FileID_ParticleDens_ColorMap, '%.16f\t%.16f\t%.16f\n', x_line(xx)/L, 0, rho(xx)*N);
      end
      fprintf(FileID_ParticleDens_ColorMap, '\n');
    end 
    

    for tt = 1:N_t
        
       %fprintf('NOTE: calculating for time=%.8f ...\n', tt*dt);   
          
       norm_diff = 9999;
       iter = 0;
          
       while (norm_diff>tol_iter) && (iter<max_iter) 
            iter = iter + 1;
            rho_old = rho;
            rho = TravellingWaves_ImpSol(rho_old,x_line,N_bgD,beta,g_bar,sigma_bar*L,tt*dt); 
            norm_diff = (1/N_x)*norm(rho_old - rho); 
       end   
          
       fprintf('SOLVING: iterative fix-point solver for time=%.8f found an %ith-step solution with norm_diff=%.16f ...\n', tt*dt, iter, norm_diff);             
       %disp(mod(tt,N_save)); %DEBUG
       
       %%% if enough timesteps is reached, printing out the results:
       %if (mod(tt,N_save)==0) || strcmp(ParticleDens_ColorMap,'yes')

            %fprintf('NOTE:printing results for time=%.8f ...\n', tt*dt);
            
            rho_print = fliplr(rho);
            rho_print = rho_print + rho - N_bgD*ones(1,N_x);
            
%            if (mod(tt,N_save)==0)
%           
%                FilenameCurrent = strcat('rho_x-LiebLiniger-GammaBG_',num2str(gamma_bg),'-dt_',num2str(dt),'-CHDsolver_TravellingWaves.t',num2str(tt*dt));
%                if exist(fullfile(DIR,FilenameCurrent),'file')
%                 fprintf('NOTE: file %s already exist; new data will be attached to its end ...\n', fullfile(DIR,FilenameCurrent));
%                 FileID_current = fopen( fullfile(DIR,FilenameCurrent) , 'at');
%                 if FileID_current==-1
%                   error('ERROR: cannot open the following file for writing, ''%s''.', fullfile(DIR,FilenameCurrent));
%                 end
%                else
%                 edit(fullfile(DIR,FilenameCurrent));
%                 FileID_current = fopen( fullfile(DIR,FilenameCurrent) , 'at');
%                 if FileID_current==-1
%                   error('ERROR: cannot open the following file for writing, ''%s''.', fullfile(DIR,FilenameCurrent));
%                 end
%                 fprintf(FileID_current,'#x/L (our unit)\t#rho*L (our unit)\n');   % printing the file header  
%                end
%
                t_array(tt+1) = tt*dt;
                for xx = 1:N_x
                   %fprintf('%.16f\t%.16f\n', x_line(xx), rho(xx));
                   %fprintf(FileID_current, '%.16f\t%.16f\n', x_line(xx)/L, rho_print(xx)*L*N_bg*(1+sqrt(2*pi)*beta*sigma_bar*erf(1/(sqrt(8.0)*sigma_bar))));
                   %full_rho(2,xx,tt+1) = x_line(xx)/L;
                   density_array(xx,tt+1) = rho_print(xx)*N;
                end

%                fclose(FileID_current);
%            
%            end
            
            if strcmp(ParticleDens_ColorMap,'yes')
              for xx = 1:N_x
                   %fprintf('%.16f\t%.16f\t%.16f\n', x_line(xx), tt*dt, rho(xx));
                   fprintf(FileID_ParticleDens_ColorMap, '%.16f\t%.16f\t%.16f\n', x_line(xx)/L, tt*dt, rho_print(xx)*N);
              end
              fprintf(FileID_ParticleDens_ColorMap,'\n');
            end  
         
     end      
       
    %end  
    
    if strcmp(ParticleDens_ColorMap,'yes')
        fclose(FileID_ParticleDens_ColorMap);
    end    

    save(fullfile(DIR,'CHD_for_tDMRG.mat'),'x_array','t_array','density_array');
    %save('CHD_for_tDMRG.mat','x_array','t_array','density_array');

    function rho = TravellingWaves_ImpSol(rho_old,x_line,N_bgD,beta,g_bar,sigma_bar,t)

      Nx = size(x_line,2);
      rho = zeros(1,Nx);

      %disp(size(rho_old)); %DEBUG
      %disp(size(rho)); %DEBUG
      %disp(size(x_line)); %DEBUG

      for xxx=1:Nx
        %%% use this one for very low gamma_bg: 
        rho(xxx) = N_bgD + 0.5*beta*N_bgD*exp( - ( x_line(xxx) - sqrt(g_bar)*(-2*sqrt(N_bgD)+3*sqrt(rho_old(xxx)))*t )^2 / (2*sigma_bar^2) );
        %%% use this one for very large gamma_bg:
        %...
      end

    end

                   
end    % the main function ends here.
