function [Mnotes, notes] = notes2matrixnotes(notes,hopsize_secs)
notes=notes(notes(:,3)>0,:); %Remove 0-pitch notes.
Hend=round(max(max(notes(end,2)),max(notes(end,1)))/hopsize_secs);
maxj=size(notes,1);
i=0;
j=1;
f0=zeros(maxj,Hend);

while (i<Hend)&&(j<=maxj)
    i=i+1;
    t = i*hopsize_secs;
    if (t>notes(j,1))&&(t<notes(j,2))
        f0(j,i)=notes(j,3);
    elseif (t>notes(j,2))
        j=j+1;
        i=i-1;
    end
end
Mnotes=f0;
notes(:,1)=notes(:,1);
notes(:,2)=notes(:,2)-notes(:,1);
notes(:,3)=notes(:,3);
end