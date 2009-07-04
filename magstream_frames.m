%%% Make movie frames from magstream
function magstream_frames(nf,ifname_root)


%%% Begin loop over frames
for frame=0:1:nf,

    %frame = 28;

    fig_name = sprintf('Movie: frame %d',frame);
    set(gcf,'name',char(fig_name));

    clf;

    %% Generate file names
    dump_file = sprintf('%s/dump%04g.dat',ifname_root,frame);
    out_file = sprintf('/home/jondata/movie/frame%04g.png',frame);

    dump_file
    out_file
    
    %% Make picture (Very important!)
    magstream(4,char(dump_file));

    set(gcf,'InvertHardcopy','off');
    set(gcf,'PaperPositionMode','auto');
        
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 750 750]);
    print(gcf,'-dpng','-r1',char(out_file));
end;