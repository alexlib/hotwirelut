function lut_calibration(calibFile,datFile,resFile,lutFile,opt)

% Usage:
% lut_calibration(calibFile,datFile,resFile,lutFile,opt)

d2r = pi/180;
% r2d = 180/pi;

% % parsing inputs: keep it here for the nexted options
% fieldNames = fieldnames(opt);
% for i = 1:length(fieldNames)
%    eval(sprintf('%s = opt.%s',fieldNames{i},fieldNames{i}));
% end


% internal variables
y_index = struct('type',[],'index',[]);
z_index = struct('type',[],'index',[]);
u_index = struct('type',[],'index',[]);

[y,z,U] = deal(zeros(opt.lenCalibFile,1));
E = zeros(opt.lenCalibFile,opt.numWires*opt.numArrays);

getCalib(calibFile);

[yg,zg] = meshgrid(min(y):opt.dy:max(y),min(z):opt.dz:max(z));

ug = interp1([u_index(:).type],linspace(1,length(u_index),opt.numelVel));

siz = [size(yg,1),size(yg,2),opt.numelVel,opt.numWires*opt.numArrays];

eg = zeros([size(yg), length(u_index), opt.numWires*opt.numArrays]);
LUT = zeros(siz);

tree = [];


if opt.tosave
    warning off
    prepareLUT();
    warning on

else
    load(lutFile)
end

if opt.testMode
    analyzeData(E');
else % dataMode

    d = readRawFld(datFile,opt.skipLength,opt.readLength);
    ss = size(d);

    [Uxhat,Uyhat,Uzhat] = deal(NaN(ss(2),opt.numArrays));
%     [yout,zout,uout] = deal(zeros(ss(2),1));

    analyzeData(d);
end

if opt.write
    writeResult(resFile);
end
% keyboard

% -----------------------------------
    function writeResult(resFile)
        %
        % keyboard
        fid = fopen(resFile,'w');
        if fid ~= -1
            fwrite(fid,[...
                Uxhat(:,1),Uyhat(:,1),Uzhat(:,1),...
                Uxhat(:,2),Uyhat(:,2),Uzhat(:,2),...
                Uxhat(:,3),Uyhat(:,3),Uzhat(:,3),...
                Uxhat(:,4),Uyhat(:,4),Uzhat(:,4),...
                Uxhat(:,5),Uyhat(:,5),Uzhat(:,5)]',...
                'float');
            fclose(fid);
        end

    end % function


% ------------------------------------------
    function prepareLUT()

        % Let's change the overall approach to array by array:
        for k = 1:opt.numArrays*opt.numWires
            for i = 1:length(u_index),
%                 eg(:,:,i,k) = gridfit(y(u_index(i).index),z(u_index(i).index),E(u_index(i).index,k),yg(1,:),zg(:,1), ...
%                     'smoothness',1, ...
%                     'interp','bilinear', ...
%                     'solver','normal', ...
%                     'regularizer','gradient', ...
%                     'extend','always', ...
%                     'tilesize',inf);
                 eg(:,:,i,k) = inpaint_nans(griddata(y(u_index(i).index),z(u_index(i).index),E(u_index(i).index,k),yg(1,:),zg(:,1)),3);

            end % i = u_index

            % keyboard
            % for m = 1:opt.numWires*opt.numArrays
            for i = 1:size(yg,1)
                for j = 1:size(yg,2)
                    LUT(i,j,1:opt.numelVel,k) = interp1([u_index(:).type],squeeze(eg(i,j,:,k)),ug','cubic');
                end
            end
        end

        tree = kdtree(reshape(LUT,[],opt.numWires*opt.numArrays));

        save(lutFile,'tree','LUT','eg');

%         warning on

    end

% --------------------------
    function getCalib(calFileName)

        z_y_U_show_plot      = 0;
        const_y_show_plot    = 0;
        const_z_show_plot    = 0;


        if nargin == 0 || isempty(calFileName)
            [y,z,U,E] = readCalib();
        else
            [y,z,U,E] = readCalib(calFileName);
        end

        % Rounding the data
        yr = round(y);
        zr = round(z);

%         Ux  =   (U).*cos(d2r*(y)).*cos(d2r*(z));
%         Uz  =   (U).*cos(d2r*(y)).*sin(d2r*(z));
%         Uy  =   (U).*sin(d2r*(y));
% fixed bug (see above): 13.08.07 
                Ux   =  U.*cos(d2r*y).*cos(d2r*z);
                Uy  =   U.*cos(d2r*z).*sin(d2r*y);
                Uz  =   U.*sin(d2r*z);


        if z_y_U_show_plot
            figure;
            plot(y,'.');hold on; plot(yr,'o');%plot all the ys
            plot(z,'r.'); plot(z,'rs'); %plot all the zs
            plot(U,'m.');%plot all the Us
            legend('y','z','U');hold off;
        end
        %

        y_index = struct('type',[],'index',[]);
        [uniqueY,junk,ind] = unique(yr);
        y_index = repmat(y_index,length(uniqueY),1);
        for i = 1:length(uniqueY)
            y_index(i).type = uniqueY(i);
            y_index(i).index = find(ind == i);
        end


        z_index = struct('type',[],'index',[]);
        [uniqueZ,junk,ind] = unique(zr);
        z_index = repmat(z_index,length(uniqueZ),1);
        for i = 1:length(uniqueZ)
            z_index(i).type = uniqueZ(i);
            z_index(i).index = find(ind == i);
        end


        %         good_monotonic_y_const = zeros(1,20);
        %         for i = 1:20
        %             for j = 1:length(y_index)
        %                 for  k = 1% :length(y_index(j).index)
        %                     good_monotonic_y_const(i) = good_monotonic_y_const(i) + ismonotonic(E(y_index(j).index(k:10:end),i));
        %                     %             plot(z(y_index(j).index(k:10:end)), E(y_index(j).index(k:10:end),i),'.--');
        %                 end
        %             end
        %         end

        %         good_monotonic_z_const = zeros(1,20);
        %         for i = 1:20
        %             for j = 1:length(z_index)
        %                 for k = 1 %:length(z_index(j).index)
        %                     good_monotonic_z_const(i) = good_monotonic_z_const(i) + ismonotonic(E(z_index(j).index(k:10:end),i));
        %                 end
        %             end
        %         end


        z_index = struct('type',[],'index',[]);
        [uniqueZ,junk,ind] = unique(zr);
        z_index = repmat(z_index,length(uniqueZ),1);
        for i = 1:length(uniqueZ)
            z_index(i).type = uniqueZ(i);
            z_index(i).index = find(ind == i);
        end

        uniqueU = mean(reshape(U,10,[]),2);
        [junk,xu] = histc(U,uniqueU+.5);
        xu = xu + 1;

        for i = 1:length(uniqueU)
            u_index(i).type = uniqueU(i);
            u_index(i).index = find(i == xu);
        end


        [n,xy] = histc(y,[y_index.type]+.5);
        [n,xz] = histc(z,[z_index.type]+.5);

        xy = xy + 1;
        xz = xz + 1;

%         %%% conical transformation
%         Uxr  =   (uniqueU(xu)).*cos(d2r*(uniqueY(xy))).*cos(d2r*(uniqueZ(xz)));
%         Uyr  =   (uniqueU(xu)).*cos(d2r*(uniqueY(xy))).*sin(d2r*(uniqueZ(xz)));
%         Uzr  =   (uniqueU(xu)).*sin(d2r*(uniqueY(xy)));
% 
%         r = sqrt(Uzr.^2 + Uyr.^2);
%         h = Uxr;
%         theta = (2*atan(r./h))*r2d;
% 
%         uniqueTheta = unique(round(theta*10)/10);
%         [junk,xtheta] = histc(theta,uniqueTheta+.5);
%         xtheta = xtheta + 1;


%         [b,ix,jx] = (unique([xtheta,xu],'rows'));


        %%%% PLOTS of Constant ys
        if const_y_show_plot
            figure
            for j = 1:20 %% for each E (in this case s(2)==20) we plot E Vs. z at const y type (in this case 9 types)
                % the plots are E vs z at const y and the graphs in the subplot are for each const velocity
                for i=1:length(y_index)
                    subplot(3,3,i); % because we have 9 types of y types
                    for k = 1:length(y_index(i).index) %every 10 places in this  array we have the same velocity
                        plot(z(y_index(i).index(k:10:end)),E(y_index(i).index(k:10:end),j),'s-');hold on;
                        title({['E','=',num2str(j),'@','\alpha','=',num2str(y_index(i).type)]})
                        xlabel('z');ylabel('E');
                    end
                    hold off; grid;
                end
            end
        end
        %close all


        %%%% PLOTS of Constant zs (the same explanation like in the above, but
        %%%% for const z types)

        if const_z_show_plot
            figure;
            for j = 1:20
                for i = 1:length(z_index)
                    subplot(4,3,i);
                    for k = 1:length(z_index(i).index)
                        plot(y(z_index(i).index(k:10:end)),E(z_index(i).index(k:10:end),j),'x-');
                        hold on;

                    end
                    hold off; grid on;
                    title({['E','=',num2str(j),'@','z','=',num2str(z_index(i).type)]});
                    xlabel('\alpha');ylabel('E');
                    %       cnt = cnt + is_monotonic(y(z_index(i).index(1:10:end)),E(j,z_index(i).index(1:10:end)));
                end
            end
        end

    end % function

% ------------------------------
    function analyzeData(d)
        if size(d,2) == opt.numWires*opt.numArrays
            d = d';
        end
        ss = size(d);
        for j = 1:ss(2)

            waitbar(j/ss(2));
%

	% disp(j)

            % looking at all wires together:
            [pntidx,pntval] = kdtree_closestpoint(tree,d(:,j)');
            [m,n,q] = ind2sub(siz,pntidx);

            % 5 x 5 neighbours, 25.07.07

            if m < 3
                m1 = 1; m2 = 5;
            elseif m > size(yg,1) - 2
                m1 = m-4; m2 = m;
            else
                m1 = m - 2;
                m2 = m + 2;
            end
            if n < 3
                n1 = 1; n2 = 5;
            elseif n > size(yg,1) - 2
                n1 = n-4; n2 = n;
            else
                n1 = n - 2;
                n2 = n + 2;
            end
            if q < 3
                q1 = 1; q2 = 5;
            elseif q > size(yg,1) - 2
                q1 = q-4; q2 = q;
            else
                q1 = q - 2;
                q2 = q + 2;
            end


            y1 = yg(m1:m2,n1:n2);
            z1 = zg(m1:m2,n1:n2);
            u1 = ug(q1:q2);

            [y2,z2,u2] = meshgrid(min(y1(:)):opt.dy2:max(y1(:)),min(z1(:)):opt.dz2:max(z1(:)),min(u1(:)):opt.du2:max(u1(:)));

            y1 = repmat(y1,[1 1 5]);
            z1 = repmat(z1,[1 1 5]);
            u1 = repmat(reshape(u1,[1 1 5]),[5 5 1]);

            %     for l = 20
            %         tmpE(:,l) = griddata3(y1,z1,u1,LUT(m-1:m+1,n-1:n+1,q-1:q+1,l),y2(:),z2(:),u2(:));
            %     end
            %
            %     for l = 1:19
            %         tmpE(:,l) = griddata3(y1,z1,u1,LUT(m-1:m+1,n-1:n+1,q-1:q+1,l),y2(:),z2(:),u2(:));
            %     end
            % l = 20
            tmpE = zeros(length(y2(:)),20);
            
            for l = 1:20
                tmpE(:,l) = myinterp3(y1,z1,u1,LUT(m1:m2,n1:n2,q1:q2,l),y2(:),z2(:),u2(:),'*cubic');
            end


            for k = 1:opt.numArrays
                wiresInArray = (1:4)+(k-1)*4; % for each array (each 4 wires)
                n2 = myNearestNeighbors([tmpE(:,wiresInArray);d(wiresInArray,j)'],1);
                % disp([k,n2(end)])
                res = [y2(n2(end)), z2(n2(end)), u2(n2(end))];
                % 13.08.07 fixed bug: in Uyhat, Uzhat.
                Uxhat(j,k)   =  (res(3)).*cos(d2r*(res(1))).*cos(d2r*(res(2)));
                Uyhat(j,k)  =   (res(3)).*cos(d2r*(res(2))).*sin(d2r*(res(1)));
                Uzhat(j,k)  =   (res(3)).*sin(d2r*(res(2)));
            end % k
        end % j

        figure
        subplot(311)
        plot(squeeze(Uxhat(:,:)));% hold on;
        subplot(312)
        plot(squeeze(Uyhat(:,:)));% hold on;
        subplot(313)
        plot(squeeze(Uzhat(:,:)));% hold on;
    end % function

% -------------------------------
    function [y,z,U,E] = readCalib(fileName)
        % [alpha,z,U,E] = READCALIB(<fileName>) reads the FLD/LOG file, using the GUI to Browse...Open
        % and returns y (horizontal angle left to right), z (vertical angle, up down), U, and E

        if nargin < 1 || isempty(fileName)
            [fileName,pathname] = uigetfile({'*.FLD;*.LOG', 'FLD or LOG file'}, 'Pick a binary .FLD or ASCII .LOG file');
        else
            [pathname,fileName,fileExt] = fileparts(fileName);
            fileName = [fileName,fileExt];
            pathname = [pathname,filesep];
        end
        % [pathstr,fileName,fileExt] = fileparts(fileName);



        if strcmp(fileName(end-3:end),'.FLD') % binary file
            d = dir([pathname,fileName]);

            fid = fopen([pathname,fileName],'r');
            if fid == -1
                error('something is wrong with the FLD file');
            end
            N = (d.bytes - 512 )/(20*2 + 3*4); % calculate the size of the calibration file
            [y,z,U] = deal(zeros(N,1));
            E = zeros(N,20);

            fseek(fid,512,'bof')
            for i = 1:N
                zzz = fread(fid, 3, 'float');
                y(i) = zzz(1);
                z(i) = zzz(2);
                U(i) = zzz(3);

                rec = fread(fid, 20, 'short');
                for k = 1:20
                    E(i,k) = rec(k);
                end
            end

        elseif strcmp(fileName(end-3:end),'.LOG')
            fid = fopen([pathname,fileName],'r');
            [y,z,U] = deal(zeros(1000,1));
            E = zeros(1000,20);
            for i = 1:16, fgetl(fid); end
            i = 0;
            while ~feof(fid)
                i = i + 1;
                tmp = str2num(fgetl(fid));
                while isempty(tmp),
                    tmp = str2num(fgetl(fid));
                end
                y(i) = tmp(1);
                z(i) = tmp(2);
                U(i) = tmp(3);
                E(i,1:20) = str2num(fgetl(fid));
                fgetl(fid);
            end
            y(i+1:end) = [];
            z(i+1:end) = [];
            U(i+1:end) = [];
            E(i+1:end,:) = [];
            fclose(fid);
        end

    end % function


% -------------------------
    function d = readRawFld(filename,skipLength,readLength)

        if nargin < 1 || isempty(filename)
            [filename,pathname] = uigetfile('*.FLD', 'Pick a RAWDATA FLD file');
        else
            [pathname,filename,fileext] = fileparts(filename);
            filename = [filename,fileext];
            pathname = [pathname,filesep];
        end

        fid = fopen([pathname,filename],'r');
        if fid == -1
            error('file does not exist %s',[pathname,filename])
        end

        nChannels = 25;

        st = fseek(fid,2048,'bof');
        if st~=0
            error('Error in positioning N');
        end

        d = fread(fid,[nChannels,readLength],'short'); % d(1:20,:) is what we need

        fclose(fid);

        d = d(1:20,:);

    end % function


end % lut1
