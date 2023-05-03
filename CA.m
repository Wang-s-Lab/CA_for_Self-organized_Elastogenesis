%%%%% CA_for_Self-organized_Elastogenesis %%%%%%%%
%%%%% Fan Xiru 20230406 %%%%%%%%%%%%%%%%%%%%%%%%%%

% EXAMPLE(with standard parameters)
% input=struct('total_Time',2400,'generate_Time',500,'n',200, ...
%              'generate_Tol',0.04,'generate_Rate',0.032, ...
%              'sigma1',0.7,'sigma2',2,'delta1',42,'delta2',38,'shape',4);
% output=struct('record',0,'radioName','ani8_150', ...
%     'saveAni',0,'strgFilename','sani401_200.mat');

function CA(input,output)
    coNei=4;%Neighbour
    %Colormap
    color_map=1;                            % background_color=[1,1,1];
    particle_color=1-[255,126,103]/255;
    cable_color=1-[7,104,159]/255;
    coSphere_color=1-[162,213,242]/255;

    computeTime =input.total_Time;
    generateTol=1.25*input.generate_Tol;
    generateRate=1.25*input.generate_Rate;
    generateTime=input.generate_Time;
    groupnolimit=input.delta1;
    groupLimit=input.delta2;
    n=input.n;

    z = zeros(n,n); 
    ani(1,computeTime).strgcells = z; 
    ani(1,computeTime).strgcoSphere = z;
    ani(1,computeTime).strgcable=z;
    cells    = z;
    cellsNew = z; 
    coSphere   =z;
    cable    =z;

    cablesNew=z;
    cableWidth=4;
    sideWidth=48;
    cable_style=input.shape;
    if cable_style==4
        cable(1:n,n-sideWidth-cableWidth:n-sideWidth)=1 ; % right ground line 
        cable(1:n,sideWidth:sideWidth+cableWidth)=1 ; % left ground line 
        cable(sideWidth:sideWidth+cableWidth, 1:n) = 1; %top line 
        cable(n-sideWidth-cableWidth:n-sideWidth, 1:n) = 1 ;%bottom line 
    elseif cable_style==3
        cable(sideWidth:sideWidth+cableWidth,1:n)=1;
        cable(n-(sideWidth+cableWidth):n-sideWidth,1:n)=1;
        for i=1:n
            if i<n-sideWidth-0.5*cableWidth
                cable(i,sideWidth+(i-1):sideWidth+cableWidth+(i-1))=1;
                cable(i,n-(sideWidth+0.5*cableWidth+(i-1)):n-(sideWidth-0.5*cableWidth+(i-1)))=1;
            else
                cable(i,i-n+sideWidth+0.5*cableWidth+1:i-n+sideWidth+1.5*cableWidth+1)=1;
                cable(i,2*n-(i+sideWidth+1.5*cableWidth):2*n-(i+sideWidth+0.5*cableWidth))=1;
            end
        end
    elseif cable_style==6
        cable(sideWidth:sideWidth+cableWidth,1.25*sideWidth:n-1.25*sideWidth)=1;
        cable(n-sideWidth-cableWidth:n-sideWidth,1.25*sideWidth:n-1.25*sideWidth)=1;
        cable(0.5*(n-cableWidth):0.5*(n+cableWidth),1:0.75*sideWidth+cableWidth)=1;
        cable(0.5*(n-cableWidth):0.5*(n+cableWidth),n-(0.75*sideWidth+cableWidth):n)=1;
        cable((n-cableWidth):n,1:0.75*sideWidth+cableWidth)=1;
        cable((n-cableWidth):n,n-(0.75*sideWidth+cableWidth):n)=1;
        for i=1:n
            if i<sideWidth+0.5*cableWidth
                cable(i,0.75*sideWidth+ceil(0.5*i)-1:0.75*sideWidth+cableWidth+ceil(0.5*i)-1)=1;
                cable(i,n-(0.75*sideWidth+cableWidth+ceil(0.5*i)-1):n-(0.75*sideWidth+ceil(0.5*i)-1))=1;
            elseif i<0.5*n
                cable(i,2.25*sideWidth+ceil(0.5*i)-1:2.25*sideWidth+cableWidth+ceil(0.5*i)-1)=1;
                cable(i,n-(2.25*sideWidth+cableWidth+ceil(0.5*i)-1):n-(2.25*sideWidth+ceil(0.5*i)-1))=1;
            elseif i<n-(sideWidth+0.5*cableWidth)
                cable(i,-0.25*sideWidth+ceil(0.5*i)-1:-0.25*sideWidth+cableWidth+ceil(0.5*i)-1)=1;
                cable(i,n-(-0.25*sideWidth+cableWidth+ceil(0.5*i)-1):n-(-0.25*sideWidth+ceil(0.5*i)-1))=1;
            else
                cable(i,1.25*sideWidth+ceil(0.5*i)-1:1.25*sideWidth+cableWidth+ceil(0.5*i)-1)=1;
                cable(i,n-(1.25*sideWidth+cableWidth+ceil(0.5*i)-1):n-(1.25*sideWidth+ceil(0.5*i)-1))=1;
            end
        end
    end
    
    sums   = z; 
    cellsums =z;
    coSums = z; 
    caSums = z;

    x = 1:n;
    y=x;
    lminus=[n,1:n-1];
    lplus=[2:n,1];

    generateParticles=1;

    if output.record==1
        v=VideoWriter(output.radioName,'MPEG-4');
        open(v);
    end
    if color_map==1
        imh=image(cat(3,1-cells*particle_color(1)-cable*cable_color(1)-coSphere*coSphere_color(1),...
                        1-cells*particle_color(2)-cable*cable_color(2)-coSphere*coSphere_color(2),...
                        1-cells*particle_color(3)-cable*cable_color(3)-coSphere*coSphere_color(3)));
    else
        imh=image(cat(3,cells,cable,coSphere));
    end
    axis equal 
    axis tight

    M=zeros(computeTime,1);

    for step=1:computeTime 
        %define neighborhood
        p=mod(step,2); %margolus neighborhood 
        %upper left cell
        xind = 1+p:2:n-2+p; 
        yind = 1+p:2:n-2+p; 
        sxind =2:1:n-1;
        syind =2:1:n-1;
        %random velocity choice 
        vary = rand(n,n)<.5 ; 
        vary1 = ~vary; 
        if step>1
            cells=ani(step-1).strgcells;
            coSphere=ani(step-1).strgcoSphere;
        end
        cellsums(sxind,syind) = cells(sxind,syind)+ cells(sxind+1,syind)+ ...
            cells(sxind,syind+1)+ cells(sxind+1,syind+1)+ cells(sxind-1,syind)+ ...
            cells(sxind,syind-1)+ cells(sxind-1,syind-1)+ cells(sxind-1,syind+1)+ ...
            + cells(sxind+1,syind-1); 
        coSphere = (cellsums>=input.sigma2) | cells>input.sigma1 | coSphere;
        if generateParticles==1
            %generate particles
            cells =cells +((rand(n,n)<generateRate)& ~cable )*generateTol;%define rate
            cells(cable==1) = 0;
            cells(coSphere==1) = 0;
            %end generation
            if step>generateTime
                    generateParticles=0;
            end
        end
        %rules of cells
        %for every particle see if it near the cable
        sums(xind,yind) = cable(xind,yind) | cable(xind+1,yind) | ...
        cable(xind,yind+1) | cable(xind+1,yind+1) ; 
        %for every particle see if it near the sphere
        coSums(sxind,syind) = coSphere(sxind+1,syind) + coSphere(sxind,syind+1) + ... 
        coSphere(sxind+1,syind+1) + coSphere(sxind-1,syind-1) + ... 
        coSphere(sxind-1,syind) + coSphere(sxind,syind-1) + ... 
        coSphere(sxind+1,syind-1) + coSphere(sxind-1,syind+1); 
        %conservation of monomers
        coSphere = ((cellsums+coSums+sums)>=4.8) | coSphere;
        cells(coSphere==1) = 0;
        %rotate the 4 cells to randomize velocity 
        cellsNew(xind,yind) = ... 
        (vary(xind,yind)  & ~sums(xind,yind)).* cells(xind+1,yind) + ... %cw 
        (vary1(xind,yind)  & ~sums(xind,yind)).* cells(xind,yind+1)+ ...
        (sums(xind,yind) .* cells(xind,yind)); %ccw 

        cellsNew(xind+1,yind) = ... 
        (vary(xind,yind)  & ~sums(xind,yind)).* cells(xind+1,yind+1)+ ... 
        (vary1(xind,yind)  & ~sums(xind,yind)).* cells(xind,yind)+...
        (sums(xind,yind) .* cells(xind+1,yind)); 

        cellsNew(xind,yind+1) = ... 
        (vary(xind,yind)  & ~sums(xind,yind)).* cells(xind,yind) + ... 
        (vary1(xind,yind)  & ~sums(xind,yind)).* cells(xind+1,yind+1)+...
        (sums(xind,yind) .* cells(xind,yind+1)); 

        cellsNew(xind+1,yind+1) = ... 
        (vary(xind,yind) & ~sums(xind,yind)).*cells(xind,yind+1) + ... 
        (vary1(xind,yind)  & ~sums(xind,yind)).* cells(xind+1,yind)+...
        (sums(xind,yind) .* cells(xind+1,yind+1)); 
        cells=cellsNew;
        ani(step).strgcells=cells;

        %rules of sphere

        %the following code divide the particles in a group into two.
        BW=coSphere;
        [~,L] = bwboundaries(BW,coNei,'noholes');
        Linitialmax=max(L(:));
        count=1;
        for i=1:max(L(:))
            [pointonum,~]=size(find(L==i));
            pointofi=find(L==i);
            if pointonum>groupnolimit
                L(pointofi(groupnolimit+1:length(pointofi)))=max(L(:))+count;
                count=count+1;
            end
        end
        BCxindex=z;
        BCyindex=z;
        movestatusL=ones(max(L(:)),1);

        xcenter=zeros(max(L(:)),1);
        ycenter=xcenter;

        for i=1:max(L(:))
            if movestatusL(i,1)==1
                point=find(L==i);%find all the point in group i
                if ~isempty(point)
                    [R,C]=ind2sub(size(coSphere),point);%find row and column of these points
                    % set the center of each group, consider Boundary problem
                    noxindex=length(find(BCxindex==i));
                    noyindex=length(find(BCyindex==i));
                    xcenter(i)=ceil((sum(R)-n*noxindex)/length(R)); %find the center of the group
                    ycenter(i)=ceil((sum(C)-n*noyindex)/length(C));
                    %throw the second part of the particles away
                    if i>Linitialmax
                        if rand>0.5
                            orientation=1;
                        else
                            orientation=-1;
                        end
                        randLen=round(rand*10);
                        xcenter(i)=xcenter(i)+randLen*orientation;
                        ycenter(i)=ycenter(i)+randLen*orientation;
                    end

                    [xcenter(i),ycenter(i)]=boundary(xcenter(i),ycenter(i),n);

                    %set the coSphere=0 and L=0 for the group first.
                    coSphere(point)=0;
                    L(point)=0;
                    %put first particle in the center
                    coSphere(xcenter(i),ycenter(i))=1;
                    L(xcenter(i),ycenter(i))=i;
                    sumxy=z;
                    %put the remain particles to the place with largest neighbors
                    if length(R)>1
                        for step1=1:length(R)-1 %reshape the group
                            sumxy(x,y)=coSphere(x,y)+coSphere(lplus,y)+coSphere(lminus,y)...
                                +coSphere(x,lplus)+coSphere(x,lminus)+coSphere(lplus,lplus)...
                                +coSphere(lplus,lminus)+coSphere(lminus,lminus)+coSphere(lminus,lplus);
                            % make the sum=0 for points who doesn't have a neighbor of group i
                            sumxy(~isneighbor(L,i))=0;
                            value = max(sumxy(:));
                            while 1
                                location=find(sumxy==value);
                                choice=find(coSphere(location)==0);
                                if ~isempty(choice)
                                    break
                                end
                                value=value-1;
                                if value==0
                                    break
                                end
                            end
                            if ~isempty(choice)
                                randa=randi(length(choice));
                                coSphere(location(choice(randa)))=1;
                                L(location(choice(randa)))=i;
                            end
                        end
                    end
                end
            end
        end
        BW=coSphere;
        [~,L] = bwboundaries(BW,coNei,'noholes');
        cableNew=z;
         %for every particle see if it near the cable
        if coNei==8
            caSums(sxind,syind) = cable(sxind,syind) + cable(sxind+1,syind) + ...
                cable(sxind,syind+1) + cable(sxind+1,syind+1) + cable(sxind-1,syind)...
                + cable(sxind,syind-1)+ cable(sxind-1,syind-1)+ cable(sxind+1,syind-1)...
                + cable(sxind-1,syind+1); %8neighbor
        elseif coNei==4
            caSums(sxind,syind) =  cable(sxind,syind)+cable(sxind+1,syind) + ...
            cable(sxind,syind+1) + cable(sxind-1,syind)+ cable(sxind,syind-1); %4neighbor
        end
        for i=1:max(L(:))
            point=find(L==i);%find all the point in group i
            length(point);
            if  length(point)>groupLimit
                movestatusL(i,1)=0;
                [R,C]=ind2sub([n,n],point);
                cableNew(R,C)=(coSphere(R,C) & (max(max(caSums(R,C)))>=3&max(max(caSums(R,C)))<=8));
            end
        end
        cablesNew=cable|cableNew|cablesNew;
        ani(step).strgcable=cablesNew;
        coSphere=coSphere & ~cableNew;
        BW=coSphere;
        [~,L] = bwboundaries(BW,coNei,'noholes');
        movestatusL=ones(max(L(:)),1);
        coSphere=randommove(coSphere,movestatusL,coNei);
        ani(step).strgcoSphere=coSphere;

        %drawnow  %need this in the loop for controls to work   
        if color_map==1
            set(imh, 'cdata', cat(3,1-cells*particle_color(1)-cablesNew*cable_color(1)-coSphere*coSphere_color(1),...
                            1-cells*particle_color(2)-cablesNew*cable_color(2)-coSphere*coSphere_color(2),...
                            1-cells*particle_color(3)-cablesNew*cable_color(3)-coSphere*coSphere_color(3)))
        else
            set(imh, 'cdata', cat(3,cells,cablesNew,coSphere))
        end
        drawnow
        if output.record==1
            M(step) = getframe;
            writeVideo(v,M(step))
        end
    end 
    if output.saveAni==1
        save (output.strgFilename,'ani','-v7.3')
    end
    if output.record==1
        close(v)
    end
end