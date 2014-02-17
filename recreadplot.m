% This program is used to help to decide the window to generate isolate filter in the gsdf code.
clear

% parameter need to be set
event='201203212215';

lalim=[20 60];
lolim=[-140 -60];
plotstanum=15;
periods=[2 20 ; 20 50];
groupv=[15 2];
prefilter=[0.5 150];
dazi=2;
sleep=0;
useramp=1;
isdist=1;


	eventinfo.name=event;
	eventinfo.gv1=groupv(1);
	eventinfo.gv2=groupv(2);
	eventinfo.t1=0;
	eventinfo.t2=1000;
	stemp=sprintf('ls %s/*.BHZ.sac > tempstalist',event);
	system(stemp);
	fpsta=fopen('tempstalist','r');
	stafile=fgetl(fpsta);
	stanum=1;
	clear stainfo stadata stadist;

	disp(' Read in sac, build up database')
	while ischar(stafile)
		% Read in sac file
		sachdr=readsac(stafile);

		if sachdr.STLA > lalim(1) && sachdr.STLA < lalim(2)...   % test whether it's in the range
				&& sachdr.STLO > lolim(1) && sachdr.STLO < lolim(2)
			stainfo(stanum).dist=sachdr.DIST;
			stainfo(stanum).stla=sachdr.STLA;
			stainfo(stanum).stlo=sachdr.STLO;
			stainfo(stanum).filename=stafile;
			stadist(stanum,1)=stanum;
			stadist(stanum,2)=sachdr.DIST;
			stadist(stanum,3)=sachdr.AZ;
			stadist(stanum,4)=sachdr.STLA;
			stadist(stanum,5)=sachdr.STLO;
			stanum=stanum+1
		end
		stafile=fgetl(fpsta);
		if ischar(stafile)
			clear sachdr
		end
	end % Loop of stations
	fclose(fpsta);

while sleep==0
	
	% Plot the world map and great circle path
	figure(4)
	clf
		evla=sachdr.EVLA; evlo=sachdr.EVLO;
		latmin=min([evla stainfo(:).stla])-10; % Bottom Latitude
		latmax=max([evla stainfo(:).stla])+10;  % Top Latitude

		longmin=min([evlo stainfo(:).stlo])-10; % Left Longitude
		longmax=max([evlo stainfo(:).stlo])+10; % Right Logitude
		
		if (longmax-longmin)>180
			temp=longmax;
			longmax=longmin+360+20;
			longmin=temp-20;
		end

		Pline=10; % Grid spacing for latidue lines
		Mline=10; % Grid spacing for longitude lines

		%define event location
		elat=evla; % Event Latitude
		elong=evlo; % Event Longitude

		hh=axesm('mapproj','aitoff',...
		   'maplatlim',[lalim(1) lalim(2)],'maplonlim',[lolim(1) lolim(2)],...
		   'MLineLocation',Mline,'PLineLocation',Pline,'Grid','on',...
		   'MeridianLabel','on','ParallelLabel','on');
		axis off, gridm on, framem on;
%		hh=worldmap('world');
		load coast
		plotm(lat,long);
		
		hold on

		plotm(elat,elong,'kx','markersize',20);
		for i=1:length(stainfo)
			plotm(stainfo(i).stla,stainfo(i).stlo,'rv')
		end
		
		%plot fine distances
	   for i=[10:10:170];
	   [azlat azlong]=scircle1(elat,elong,i);
	   geoshow(azlat,azlong,'color','k','LineStyle','--','LineWidth',.5);
	   end

		[epidist baz]=distance(mean([stainfo(:).stla]),mean([stainfo(:).stlo]),evla,evlo);
		stemp=sprintf('%s \n Dist=%f, Baz=%f Mag:%3f',event,epidist,baz,sachdr.MAG);
		title(stemp);

		[mla mlo]=inputm(2);

		[stamla stamlo]=gcwaypts(mla(1),mlo(1),mla(2),mlo(2),plotstanum);

%		mazi=azimuth(elat,elong,mla,mlo);
	% End of plot the world map
	

	% sort the station by epicenter distance
	stadist_sort=sortrows(stadist,2);

	% select some stations to plot
	stanum=length(stadist);
%	ddist=ceil((stadist_sort(end,2)-stadist_sort(1,2))/plotstanum);
%	dist=stadist_sort(1,2):ddist:stadist_sort(end,2);  % dist is the array that stations being plotted.

	% Find and store the stations being plotted in the structure stadata
	disp('Finding the station to pick\n')
	maxamp=0;
	for i=1:plotstanum
		disp(i/plotstanum);
		[mdist stai] = min(abs(distance(stamla(i),stamlo(i),stadist(:,4),stadist(:,5))));
		sac=readsac(stainfo(stai).filename);
		stadata(i).filename=stainfo(stai).filename;
%		stadata(i).data=sac.DATA1';
		f1=1/prefilter(2);
		f2=1/prefilter(1);
		fN=1/2/sac.DELTA;
		[b,a]=butter(2,[f1/fN, f2/fN]);
		stadata(i).data=filter(b,a,sac.DATA1');
		if maxamp < max(abs((stadata(i).data)))
			maxamp = max(abs((stadata(i).data)));
		end
		stadata(i).dist=sac.DIST;
		stadata(i).az=sac.AZ;
		stadata(i).stla=sac.STLA;
		stadata(i).stlo=sac.STLO;
		stadata(i).t=sac.B:sac.DELTA:sac.B+(sac.NPTS-1)*sac.DELTA;
	end
	for i=1:plotstanum
		stadata(i).data=stadata(i).data./maxamp; % Normalize it
	end
	dist(1)=min([stadata(:).dist]);
	dist(2)=max([stadata(:).dist]);
	az(1)=min([stadata(:).az]);
	az(2)=max([stadata(:).az]);
	ddist=(dist(2)-dist(1))/plotstanum;
	daz=(az(2)-az(1))/plotstanum;
	eventinfo.gv1=groupv(1);
	eventinfo.gv2=groupv(2);
	eventinfo.t1=0;
	eventinfo.t2=1000;	

	% filter the waveforms
	disp('Filter the waveforms');
	fN=1/2/sac.DELTA;
	maxampf=zeros(length(periods),1);
	for i=1:length(stadata)
		disp(i/plotstanum);
		for j=1:size(periods,1)
			f1=1/periods(j,2);
			f2=1/periods(j,1);
			[b,a]=butter(2,[f1/fN, f2/fN]);
			stadata(i).fdata(j,:)=filter(b,a,stadata(i).data);
			if maxampf(j) < max(abs(stadata(i).fdata(j,:)))
				maxampf(j) = max(abs(stadata(i).fdata(j,:)));
			end
%			stadata(i).fdata(j,:)=stadata(i).fdata(j,:)./max(abs(stadata(i).fdata(j,:)));
		end
	end

	% Normalize filtered data
	for i=1:length(stadata)
		for j=1:size(periods,1)
			stadata(i).fdata(j,:)=stadata(i).fdata(j,:)./maxampf(j);
		end
	end

	bgx=dist(1)/groupv(1)-1000;
	endx=dist(end)/groupv(2)+1000;
	% Plot the original data
	figure(1)
	clf
		hold on
		if isdist
			amp=(dist(end)-dist(1))/plotstanum*useramp;
		else
			amp=daz*useramp;
		end
		for i=1:length(stadata)
			if isdist
				plot(stadata(i).t, stadata(i).data*amp+stadata(i).dist);
			else
				plot(stadata(i).t, stadata(i).data*amp+stadata(i).az);
			end
		end
		if isdist
			plot([dist(1)/groupv(1) dist(end)/groupv(1)],[dist(1) dist(end)],'r')
			plot([dist(1)/groupv(2) dist(end)/groupv(2)],[dist(1) dist(end)],'r')
			ylim([dist(1)-ddist dist(end)+ddist])
			ylabel('Dist /km')
		else
			ylabel('Azimuth /km')

		end
		xlim([bgx endx])
		% Plot filtered data
		%
	figure(2)
	clf
	gv1=groupv(1);gv2=groupv(2);
	t1=0;t2=0;
	for i=1:length(periods)
		subplot(1,length(periods),i)
			hold on
			if isdist
				amp=(dist(end)-dist(1))/plotstanum*useramp;
			else
				amp=daz*useramp;
			end

			for j=1:length(stadata)
				if isdist
					plot(stadata(j).t, stadata(j).fdata(i,:)*amp+stadata(j).dist);
				else
					plot(stadata(j).t, stadata(j).fdata(i,:)*amp+stadata(j).az);
				end
			end
			if isdist
				plot([dist(1)/gv1+t1 dist(end)/gv1+t1],[dist(1) dist(end)],'r')
				plot([dist(1)/gv2+t2 dist(end)/gv2+t2],[dist(1) dist(end)],'r')
				ylabel('Dist /km')
			else
				ylabel('Azimuth /km')
			end

			if isdist
				ylim([dist(1)-ddist dist(end)+ddist])
			end
			xlim([bgx endx])
			title(sprintf('%d - %d s',periods(i,1),periods(i,2)));
	end

	% Plot the world map and great circle path
	figure(3)
	clf
		evla=sac.EVLA; evlo=sac.EVLO;
		latmin=min([evla stadata(:).stla])-10; % Bottom Latitude
		latmax=max([evla stadata(:).stla])+10;  % Top Latitude

		longmin=min([evlo stadata(:).stlo])-10; % Left Longitude
		longmax=max([evlo stadata(:).stlo])+10; % Right Logitude
		
		if (longmax-longmin)>180
			temp=longmax;
			longmax=longmin+360+20;
			longmin=temp-20;
		end

		Pline=10; % Grid spacing for latidue lines
		Mline=10; % Grid spacing for longitude lines

		%define event location
		elat=evla; % Event Latitude
		elong=evlo; % Event Longitude

		hh=axesm('mapproj','aitoff',...
		   'maplatlim',[latmin latmax],'maplonlim',[longmin longmax],...
		   'MLineLocation',Mline,'PLineLocation',Pline,'Grid','on',...
		   'MeridianLabel','on','ParallelLabel','on');
		axis off, gridm on, framem on;
%		hh=worldmap('world');
		load coast
		plotm(lat,long);
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(hh, states, 'FaceColor', [0.5 0.5 1])
		
		hold on

		for i=1:length(stadata)
		   [la, lo]=gcwaypts(elat,elong,stadata(i).stla,stadata(i).stlo,30);
		   hh=geoshow(la,lo,'displaytype','line','color','r');
		   set(hh,'LineWidth',1)
		end
		
		[epidist baz]=distance(mean([stadata(:).stla]),mean([stadata(:).stlo]),evla,evlo);
		stemp=sprintf('%s \n Dist=%f, Baz=%f Mag:%3f',event,epidist,baz,sachdr.MAG);
		title(stemp);
	% End of plot the world map

	% Interact with user
	while 1
		disp('What do you want to do?')
		disp('1: Change time axis')
		disp('2: Change mark')
		disp('3: Change between dist and azi')
		disp('4: reset the window')
		disp('5: reselect the station');
		disp('6: change the amplitude');
		disp('7: quit')
		resp=input('','s');
		if resp=='7'
			disp('Ok, quit')
			sleep=1;
			break;
		end
		if resp=='6'
			useramp=input('','s');
			useramp=str2num(useramp);
		end
		if resp=='5'
			disp('Ok, reselect')
			sleep=0;
			break;
		end
		if resp=='4'
			disp('Ok, reset the window')
			eventinfo.gv1=groupv(1);
			eventinfo.gv2=groupv(2);
			eventinfo.t1=0;
			eventinfo.t2=1000;
			dist(1)=min([stadata(:).dist]);
			dist(2)=max([stadata(:).dist]);
			useramp=1;
		end
		if resp=='1'
				figure(1)
				coor=ginput(2);
				bgx=min(coor(:,1));
				endx=max(coor(:,1));
%				A=[coor(1,2) 1;coor(2,2) 1];
%				B=[coor(1,1);coor(2,1)];
%				x=A\B;
%				eventinfo.gv1=1/x(1);
%				eventinfo.t1=x(2);
		end
		if resp=='2'
			if isdist
				figure(1)
				coor=ginput(2);
				A=[coor(1,2) 1;coor(2,2) 1];
				B=[coor(1,1);coor(2,1)];
				x=A\B;
				eventinfo.gv2=1/x(1);
				eventinfo.t2=x(2);
			else
				disp('Please change to dist plot first');
			end
		end
		if resp=='3'
			isdist=~isdist;
		end

		disp([eventinfo.gv1 eventinfo.t1]);
		disp([eventinfo.gv2 eventinfo.t2]);

		% Now, we need to flash the waveform plot by new window
		figure(1)
		clf
			hold on
			if isdist
				amp=(dist(end)-dist(1))/plotstanum*useramp;
			else
				amp=daz*useramp;
			end
			for i=1:length(stadata)
				if isdist
					plot(stadata(i).t, stadata(i).data*amp+stadata(i).dist);
				else
					plot(stadata(i).t, stadata(i).data*amp+stadata(i).az);
				end
			end
			gv1=eventinfo.gv1;
			gv2=eventinfo.gv2;
			t1=eventinfo.t1;
			t2=eventinfo.t2;
			if isdist
				plot([dist(1)/gv1+t1 dist(end)/gv1+t1],[dist(1) dist(end)],'r')
				plot([dist(1)/gv2+t2 dist(end)/gv2+t2],[dist(1) dist(end)],'r')
				ylabel('Dist /km')
			else
				ylabel('Azimuth /km')
			end
			if isdist
				ylim([dist(1)-ddist dist(end)+ddist])
			end
			xlim([bgx endx]);

		% Plot filtered data
		figure(2)
		clf
		for i=1:length(periods)
			subplot(1,length(periods),i)
				hold on
				if isdist
					amp=(dist(end)-dist(1))/plotstanum*useramp;
				else
					amp=daz*useramp;
				end

				for j=1:length(stadata)
					if isdist
						plot(stadata(j).t, stadata(j).fdata(i,:)*amp+stadata(j).dist);
					else
						plot(stadata(j).t, stadata(j).fdata(i,:)*amp+stadata(j).az);
					end
				end
				if isdist
					plot([dist(1)/gv1+t1 dist(end)/gv1+t1],[dist(1) dist(end)],'r')
					plot([dist(1)/gv2+t2 dist(end)/gv2+t2],[dist(1) dist(end)],'r')
					ylim([dist(1)-ddist dist(end)+ddist])
					ylabel('Dist /km')
				else
					ylabel('Azimuth /km')
				end
				xlim([bgx endx]);
				xlabel('time /s')
				title(sprintf('%d - %d s',periods(i,1),periods(i,2)));
		end

	end % End of input loop

end

