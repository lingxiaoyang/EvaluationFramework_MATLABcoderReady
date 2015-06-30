function [Mnotes, notes] = midi2matrixnotes(midifile,hopsize_secs)
midi=readmidi(midifile);
M= midiInfo(midi);
M=M(M(:,3)>0,:); %Remove 0-pitch notes.
Hend=round(max(max(M(end,6)),max(M(end,5)))/hopsize_secs);
maxj=size(M,1);
i=0;
j=1;
f0=zeros(maxj,Hend);

while (i<Hend)&&(j<=maxj)
    i=i+1;
    t = i*hopsize_secs;
    if (t>M(j,5))&&(t<M(j,6))
        f0(j,i)=M(j,3);
    elseif (t>M(j,6))
        j=j+1;
        i=i-1;
    end
end
Mnotes=f0;
notes(:,1)=M(:,5);
notes(:,2)=M(:,6)-M(:,5);
notes(:,3)=M(:,3);
end