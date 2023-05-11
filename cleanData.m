clear
fclose('all');

cleanYr = true;
station = {'USM00074494','USM00072501','USM00072694'};

if cleanYr
    for i=1:numel(station)

        writeFlag = false;

        fidRead = fopen(['Data\' station{i} '-data.txt']);
        fidWrite = fopen(['Data\' station{i} '-data-beg2021.txt'],'w');
        while ~feof(fidRead)
            tmpline = fgetl(fidRead);
            if (tmpline(1)=='#')
                if (tmpline(25:26)=='00') % only log data collected at midnight
                    writeFlag = true;
                else
                    writeFlag = false;
                end
            end
            if writeFlag
                if (tmpline(11)=='-') % check if PRESS exists
                    continue
                elseif (tmpline(17)=='-') % check if GPH exists
                    continue
                elseif (tmpline(23)=='-') % check if TEMP exists
                    continue
                elseif (tmpline(41)=='-') % check if WDIR exists
                    continue
                elseif (tmpline(47)=='-') % check if WSPD exists
                    continue
                end
                fprintf(fidWrite,'%s\r\n',tmpline);
            end
        end
        fclose(fidWrite);
        fclose(fidRead);

    end
end
