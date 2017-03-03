% Suppose in a plant of size 30m by 60m, many sensors are installed on the 
% ceiling to monitor the environmental temperature. The sensors communicate 
% via wireless data transmission. Suppose the transmission radius is 3m. 
% Assume that there is no central unit to collect and process the data.
% All computations and transmissions are carried over the sensors 
% in a decentralized way. 

  close all;
  clear;
  clc;

  r = 3;
  l = 60;
  w = 30;
  

% Sensors placed as a 2D grid

  N = 0;      % sensor number
  row = 0;    % tmp variables for placing sensors
  column = 0; % 

  dy = r;     % distance between the sensors in both directions
  dx = r;           
  y_margin = r/sqrt(2);    % distance from the wall
  x_margin = r/sqrt(2);     

 % figure('Position',[600 300 600 500]);
%  rectangle('Position',[0,0,l,w]);    % walls of the plant
%  axis([-5 65 -5 35]);                       
 % xlabel('Length [m]'); ylabel('Width [m]');

  while (column-1)*dx+2*x_margin < l
    figure(99)
      hold on;

    while (row-1)*dy+2*y_margin < w
      x = x_margin+column*dx;
      y = row*dy+y_margin;
      N = N+1;
      coordinates(N,:) = [x y]; 
      plot(x,y,'o');
      %circle(x,y,3);
      row=row+1;
    end

    row=0;
    column=column+1;
  end
  axis([0 60 0 30]);

 % title(sprintf('Plant: (%i m x %i m) with %i sensors', l, w, N));

  
% For all N nodes find all nodes that are in range (adjacent)
% in case of 2D grid 4 nodes will have 2 connected, boundary nodes - 
% 3 connected, and all others - 4 nodes connected
% but it is handy to know their locations for organizing information fusion

% build an adjacency matrix of the sensor grid
% if two nodes are different and the distance between them
% is smaller than the range - then these nodes are adjacent

adjacent = zeros(N,N);
  for i=1:N
    for j=1:N
      if (i ~= j) && (norm(coordinates(i,:) - coordinates(j,:),2) <= r) 
        adjacent(i,j) = 1;
      end
    end
  end
  
N_adj = sum(adjacent);

    
    
 % Averaging Problem Solution using ADMM  
% solve the problem described on slide 61
% by realizing an update algorithm on slide 62


% z and v have dimensions Nx4 because we need one fusion center per each
% edge and for 2D grid node has maximum 4 links connected

% an algorithm is as follows:
% 
a_i = 20+3.*randn(N,1);% need to assign random values (sensor measurements)
x_i = zeros(N,1); % most recent estimate of the mean
z = zeros(N,4);
v = zeros(N,4);
x_n = zeros(N,4); %x values for neighbors
n_ind = zeros(N,4); % indices for neighbors
rho = 1;
% all info about nodes and their measurements
% each node have access to only the info on its corresponding row
% data = |coordinates (Nx2)|a_i (Nx1)|x_i (Nx1)|z (Nx4)|v(Nx4)|x_n(Nx4)|
%        |n_ind (Nx4)|
%        |====== 1:2 ======|== 3:3 ==|== 4:4 ==|= 5:8 =| 9:12 | 13:16 =|
%        |== 17:20 ==|
% each node needs to have info of adjacent xi's as well!
foo = 1:N;
for i = 1:N
    %x_n(i,:) = horzcat(x_i(adjacent(:,i)==1)',zeros(1,4-length(x_i(adjacent(:,i)==1))));
    n_ind(i,:) = horzcat(foo(adjacent(:,i)==1),zeros(1,4-length(foo(adjacent(:,i)==1))));
end
data = horzcat(coordinates, a_i, x_i, z, v, x_n, n_ind);
run_type = input(['Choose run type:\n0 - Synchronous, fixed iterations\n1 - Synchronous, fixed error' ...
    '\n2 - Asynchronous, fixed iteration\n3 - Asynchronous, fixed error\n']);
switch(run_type)
    case(0)
        %% Synchronous - fixed iteration
        p = 0:0.5:4.5; % prob. of package loss = qfunc(p)
        pr_loss = p(input('Probability of error:\n1 (most, p=0.5) - 10 (least, p = 4e-6)\n'));
        %[1-(maxErr+1)] transmission errors in random sensors - uniform dist
        maxErr = 10:5:55;
        faulty_sens_N = maxErr(input('Number of faulty sensors (1-10)\n'));
        total_trans = 0;
        loss_type = input('Loss type (0 - xnew(k+1)=xnew(k); 1 - xnew(k+1)=0)\n');
        max_iter = 100;
        avgs = zeros(N,max_iter);
        err = zeros(1,max_iter);
        faults=[];
        for i=1:max_iter
            xnew = (data(:,3)+sum(rho*data(:,5:8)-data(:,9:12),2))./(1+rho*N_adj)';
            if (randn(1,1)>pr_loss) 
                Loss = 1;
                else Loss = 0;
            end
            if(Loss)
                faults = round(1+(N-1).*rand(ceil(1+faulty_sens_N*rand(1,1)),1));
                if(loss_type) 
                    xnew(faults) = 0;
                else
                    xnew(faults) = data(faults,3);
                end
            end
            data(:,4) = xnew;
            znew = zeros(N,4);
            vnew = zeros(N,4);
            for j = 1:N
            neighs = nonzeros(data(j,17:20));
            v_ji = zeros(1,4);
            xj = zeros(1,4);
            k=1;
            %find v_ji
            for k=1:N_adj(j)
                ind = find(data(neighs(k),17:20)==j);
            % each node transmits v_ji to its neighbors (Nxk trans)
                v_ji(k) = data(neighs(k),ind+8);
            % each node transmits x_j to its neigbors (Nxk trans)
                xj(k) = xnew(neighs(k)); 
            end
            total_trans = total_trans + 2*N*k;
            znew(j,:)= 0.5*(rho^-1*(data(j,9:12)+...    %v_ij
                v_ji)+...                               %v_ji
                repmat(xnew(j),1,4)...                  %x_i
                +xj);                                   %x_j
            vnew(j,:) = data(j,9:12)+rho*(repmat(xnew(j),1,4)-znew(j,:));
            end
            if(Loss)
                if(loss_type)
                    znew(faults,:) = zeros(length(faults),4);
                    vnew(faults,:) = zeros(length(faults),4); 
                else
                    znew(faults,:) = data(faults,5:8);
                    vnew(faults,:) = data(faults,5:8);
                end
            end

            total_trans = total_trans + N*k;
            data(:,5:8) = znew;
            data(:,9:12) = vnew;
            avgs(:,i) = data(:,4);
            err(i) = norm(mean(data(:,3))-mean(avgs(:,i)),1)/mean(data(:,3))*100;
            xnew = zeros(N,1);
            if(i==30)
                data(:,3) = 30+3.*randn(N,1);
            end
        end
        figure(2)
        plot(err)
        ylabel('Percent error');
        xlabel('Number of iterations');
        figure(3)
        % plot some of the averages
        plot(avgs(randperm(N,30),:)'); % taking 30 random nodes
        hold on;
        %title('Average values stored in all N sensors (syncronous, no faults)');
        xlabel('Number of iterations');
        ylabel('Average temperature values');
    case(1)
        %% Synchronous - fixed error
        p = linspace(1,4);
        %[1-(maxErr+1)] transmission errors in random sensors - uniform dist
        maxErr = 10:5:55;
        faulty_sens_N = maxErr(input('Number of faulty sensors (1-10)\n'));
        eps = input('Error to stop at\n');
        loss_type = input('Loss type (0 - xnew(k+1)=xnew(k); 1 - xnew(k+1)=0)\n');
        iter = zeros(1,length(p));
        reverseStr = '';
        for i = 1:length(p)
            a_i = 20+3.*randn(N,1);
            data(:,3) = a_i;
            data(:,4:12) = zeros(N,9);
            err = 100;
            while(err>eps)
                xnew = (data(:,3)+sum(rho*data(:,5:8)-data(:,9:12),2))./(1+rho*N_adj)';
                if (randn(1,1)>p(i)) 
                    Loss = 1;
                    else Loss = 0;
                end
                if(Loss)
                    if(loss_type)
                        faults = round(1+(N-1).*rand(ceil(1+faulty_sens_N*rand(1,1)),1));
                        xnew(faults) = 0;
                    else
                        faults = round(1+(N-1).*rand(ceil(1+faulty_sens_N*rand(1,1)),1));
                        xnew(faults) = data(faults,3);
                    end
                end
                data(:,4) = xnew;
                znew = zeros(N,4);
                vnew = zeros(N,4);
                for j = 1:N
                    neighs = nonzeros(data(j,17:20));
                    v_ji = zeros(1,4);
                    xj = zeros(1,4);
                    k=1;
                    %find v_ji
                    for k=1:N_adj(j)
                        ind = find(data(neighs(k),17:20)==j);
                        % each node transmits v_ji to its neighbors (Nxk trans)
                        v_ji(k) = data(neighs(k),ind+8);
                        % each node transmits x_j to its neigbors (Nxk trans)
                        xj(k) = xnew(neighs(k)); 
                    end
                        znew(j,:)= 0.5*(rho^-1*(data(j,9:12)+...    %v_ij
                        v_ji)+...                               %v_ji
                        repmat(xnew(j),1,4)...                  %x_i
                        +xj);                                   %x_j
                    vnew(j,:) = data(j,9:12)+rho*(repmat(xnew(j),1,4)-znew(j,:));
                end
                if(Loss)
                    if(loss_type)
                        znew(faults,:) = zeros(length(faults),4); 
                        vnew(faults,:) = zeros(length(faults),4);                                                       
                    else
                        znew(faults,:) = data(faults,5:8);
                        vnew(faults,:) = data(faults,5:8);
                    end
                end                     
                data(:,5:8) = znew;
                data(:,9:12) = vnew;
                err = norm(mean(data(:,3))-mean(data(:,4)),1)/mean(data(:,3))*100;
                xnew = zeros(N,1);
                iter(i) = iter(i) + 1;
                msg = sprintf('%d/%d', i, length(p));
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
        fprintf('\n');
        figure(4)
        plot(qfunc(p),iter);
        title('Packet loss probability vs number of iterations required');
        ylabel('Number of iterations');
        xlabel('Packet loss probability (per iteration)');
    case(2)
        %% Asynchronous - fixed iterations
        % One part of the network calculates and updates its values first, then the
        % second part does the same.
        div = input('Network division factor [0-1]\n');
        N1 = round(N*div);
        N2 = N-N1;
        p = 0:0.5:4.5; % prob. of package loss = qfunc(p)
        pr_loss = p(input('Probability of error:\n1 (most, p=0.5) - 10 (least, p = 4e-6)\n'));
        %[1-(maxErr+1)] transmission errors in random sensors - uniform dist
        maxErr = 10:5:55;
        faulty_sens_N = maxErr(input('Number of faulty sensors (1-10)\n'));
        total_trans = 0;
        loss_type = input('Loss type (0 - xnew(k+1)=xnew(k); 1 - xnew(k+1)=0)\n');
        max_iter = 100;
        avgs = zeros(N,max_iter);
        err = zeros(1,max_iter);
        for i=1:2*max_iter
            if (randn(1,1)>pr_loss) 
                Loss = 1;
                else Loss = 0;
            end
            % First part
            if(mod(i,2)==0)
                xnew1 = (data(1:N1,3)+sum(rho*data(1:N1,5:8)-data(1:N1,9:12),2))./(1+rho*N_adj(1:N1))';
                if(Loss)
                    faults = round(1+(N1-1).*rand(ceil(1+faulty_sens_N*rand(1,1)),1));
                    if(loss_type)                        
                        xnew1(faults) = 0;
                    else
                        xnew1(faults) = data(faults,3);
                    end
                end
                data(1:N1,4) = xnew1;
                znew1 = zeros(N1,4);
                vnew1 = zeros(N1,4);
                for j = 1:N1
                    neighs = nonzeros(data(j,17:20));
                    v_ji1 = zeros(1,4);
                    xj1 = zeros(1,4);
                    k=1;
                    %find v_ji
                    for k=1:N_adj(j)
                        ind = find(data(neighs(k),17:20)==j);
                        %each node transmits v_ji to its corresponding neighbors (Nxk
                        %transmissions)
                        v_ji1(k) = data(neighs(k),ind+8);
                        %each node transmits x_j values to its neigbors (Nxk
                        %transmissions)
                        xj1(k) = data(neighs(k),4); 
                    end
                    znew1(j,:)= 0.5*(rho^-1*(data(j,9:12)+...   %v_ij
                        v_ji1)+...                              %v_ji
                        repmat(xnew1(j),1,4)...                 %x_i
                        +xj1);                                  %x_j
                    vnew1(j,:) = data(j,9:12)+rho*(repmat(xnew1(j),1,4)-znew1(j,:));
                end
                if(Loss)
                    if(loss_type)
                        znew1(faults,:) = zeros(length(faults),4); 
                        vnew1(faults,:) = zeros(length(faults),4);                                                       
                    else
                        znew1(faults,:) = data(faults,5:8);
                        vnew1(faults,:) = data(faults,5:8);
                    end
                end  
                data(1:N1,5:8) = znew1;
                data(1:N1,9:12) = vnew1;
                xnew1 = zeros(N1,1);
            else
            % Second part
                xnew2 = (data(N1+1:N,3)+sum(rho*data(N1+1:N,5:8)-data(N1+1:N,9:12),2))./(1+rho*N_adj(N1+1:N))';
                if(Loss)
                    faults = round(1+(N2-1).*rand(ceil(1+faulty_sens_N*rand(1,1)),1));
                    if(loss_type)
                        xnew2(faults) = 0;
                    else
                        xnew2(faults) = data(faults,3);
                    end
                end
                data(N1+1:N,4) = xnew2;
                znew2 = zeros(N2,4);
                vnew2 = zeros(N2,4);
                for j = N1+1:N
                    neighs = nonzeros(data(j,17:20));
                    v_ji2 = zeros(1,4);
                    xj2 = zeros(1,4);
                    k=1;
                    %find v_ji
                    for k=1:N_adj(j)
                        ind = find(data(neighs(k),17:20)==j);
                        %each node transmits v_ji to its corresponding neighbors (Nxk
                        %transmissions)
                        v_ji2(k) = data(neighs(k),ind+8);
                        %each node transmits x_j values to its neigbors (Nxk
                        %transmissions)
                        xj2(k) = data(neighs(k),4); 
                    end
                    %total_trans = total_trans + 2*N2*k;
                    znew2(j-N1,:)= 0.5*(rho^-1*(data(j,9:12)+...    %v_ij
                        v_ji2)+...                                  %v_ji
                        repmat(data(j,4),1,4)...                    %x_i
                        +xj2);                                      %x_j
                    vnew2(j-N1,:) = data(j,9:12)+rho*(repmat(data(j,4),1,4)-znew2(j-N1,:));
                end
                if(Loss)
                    if(loss_type)
                        znew2(faults,:) = zeros(length(faults),4); 
                        vnew2(faults,:) = zeros(length(faults),4);                                                       
                    else
                        znew2(faults,:) = data(faults,5:8);
                        vnew2(faults,:) = data(faults,5:8);
                    end
                end  
                data(N1+1:N,5:8) = znew2;
                data(N1+1:N,9:12) = vnew2;
                xnew2 = zeros(N2,1);
            end
            avgs(:,i) = data(:,4);
            err(i) = norm(mean(data(:,3))-mean(avgs(:,i)),1)/mean(data(:,3))*100;
        end
        figure(5)
        plot(err)
        %title('Percent error from true average');
        ylabel('Percentage');
        xlabel('Number of iterations');
        figure(6)
        %plot some of the averages
        plot(avgs(randperm(N,30),:)'); % taking 30 random nodes
        %title('Average values stored in sensors 30-50');
        xlabel('Number of iterations');
        ylabel('Average temperature values');
        %data(:,3:12) = zeros(N,10);
    case(3)
        %% Asynchronous - fixed error
        div = input('Network division factor [0-1]\n');
        N1 = round(N*div);
        N2 = N-N1;
        p = linspace(1,4);
        %[1-(maxErr+1)] transmission errors in random sensors - uniform dist
        maxErr = 10:5:55;
        faulty_sens_N = maxErr(input('Number of faulty sensors (1-10)\n'));
        eps = input('Error to stop at\n');
        loss_type = input('Loss type (0 - xnew(k+1)=xnew(k); 1 - xnew(k+1)=0)\n');
        iter = zeros(1,length(p));
        reverseStr = '';
        for i = 1:length(p)
            a_i = 20+3.*randn(N,1);
            data(:,3) = a_i;
            data(:,4:12) = zeros(N,9);
            err = 100;
            while(err>eps)
                if (randn(1,1)>p(i)) 
                    Loss = 1;
                    else Loss = 0;
                end
                if(mod(iter(i),2)==0)
                    xnew1 = (data(1:N1,3)+sum(rho*data(1:N1,5:8)-data(1:N1,9:12),2))./(1+rho*N_adj(1:N1))';
                    if(Loss)
                        faults = round(1+(N1-1).*rand(ceil(1+faulty_sens_N*rand(1,1)),1));
                        if(loss_type)
                            xnew1(faults) = 0;
                        else
                            xnew1(faults) = data(faults,3);
                        end
                    end
                    data(1:N1,4) = xnew1;
                    znew1 = zeros(N1,4);
                    vnew1 = zeros(N1,4);
                    for j = 1:N1
                        neighs = nonzeros(data(j,17:20));
                        v_ji1 = zeros(1,4);
                        xj1 = zeros(1,4);
                        k=1;
                        %find v_ji
                        for k=1:N_adj(j)
                            ind = find(data(neighs(k),17:20)==j);
                            %each node transmits v_ji to its corresponding neighbors (Nxk
                            %transmissions)
                            v_ji1(k) = data(neighs(k),ind+8);
                            %each node transmits x_j values to its neigbors (Nxk
                            %transmissions)
                            xj1(k) = data(neighs(k),4); 
                        end
                        znew1(j,:)= 0.5*(rho^-1*(data(j,9:12)+...   %v_ij
                            v_ji1)+...                              %v_ji
                            repmat(xnew1(j),1,4)...                 %x_i
                            +xj1);                                  %x_j
                        vnew1(j,:) = data(j,9:12)+rho*(repmat(xnew1(j),1,4)-znew1(j,:));
                    end
                    if(Loss)
                        if(loss_type)
                            znew1(faults,:) = zeros(length(faults),4); 
                            vnew1(faults,:) = zeros(length(faults),4);                                                       
                        else
                            znew1(faults,:) = data(faults,5:8);
                            vnew1(faults,:) = data(faults,5:8);
                        end
                    end  
                    data(1:N1,5:8) = znew1;
                    data(1:N1,9:12) = vnew1;
                    xnew1 = zeros(N1,1);
                else
                    % Second part
                    xnew2 = (data(N1+1:N,3)+sum(rho*data(N1+1:N,5:8)-data(N1+1:N,9:12),2))./(1+rho*N_adj(N1+1:N))';
                    if(Loss)
                        faults = round(1+(N2-1).*rand(ceil(1+faulty_sens_N*rand(1,1)),1));
                        if(loss_type)
                            xnew2(faults) = 0;
                        else
                            xnew2(faults) = data(faults,3);
                        end
                    end
                    data(N1+1:N,4) = xnew2;
                    znew2 = zeros(N2,4);
                    vnew2 = zeros(N2,4);
                    for j = N1+1:N
                        neighs = nonzeros(data(j,17:20));
                        v_ji2 = zeros(1,4);
                        xj2 = zeros(1,4);
                        k=1;
                        %find v_ji
                        for k=1:N_adj(j)
                            ind = find(data(neighs(k),17:20)==j);
                            %each node transmits v_ji to its corresponding neighbors (Nxk
                            %transmissions)
                            v_ji2(k) = data(neighs(k),ind+8);
                            %each node transmits x_j values to its neigbors (Nxk
                            %transmissions)
                            xj2(k) = data(neighs(k),4); 
                        end
                        %total_trans = total_trans + 2*N2*k;
                        znew2(j-N1,:)= 0.5*(rho^-1*(data(j,9:12)+...    %v_ij
                            v_ji2)+...                                  %v_ji
                            repmat(data(j,4),1,4)...                    %x_i
                            +xj2);                                      %x_j
                        vnew2(j-N1,:) = data(j,9:12)+rho*(repmat(data(j,4),1,4)-znew2(j-N1,:));
                    end
                    if(Loss)
                        if(loss_type)
                            znew2(faults,:) = zeros(length(faults),4); 
                            vnew2(faults,:) = zeros(length(faults),4);                                                       
                        else
                            znew2(faults,:) = data(faults,5:8);
                            vnew2(faults,:) = data(faults,5:8);
                        end
                    end  
                    data(N1+1:N,5:8) = znew2;
                    data(N1+1:N,9:12) = vnew2;
                    xnew2 = zeros(N2,1);
                end
                err = norm(mean(data(:,3))-mean(data(:,4)),1)/mean(data(:,3))*100;
                iter(i) = iter(i) + 1;
                msg = sprintf('%d/%d', i, length(p));
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
        fprintf('\n');
        figure(4)
        plot(qfunc(p),iter);
        title('Packet loss probability vs number of iterations required');
        ylabel('Number of iterations');
        xlabel('Packet loss probability (per iteration)');
        
end