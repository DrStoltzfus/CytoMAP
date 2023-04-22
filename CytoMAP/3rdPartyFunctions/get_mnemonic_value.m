
function mneval = get_mnemonic_value(mnemonic_name,fcsheader,mnemonic_separator)

    if strcmp(mnemonic_separator,'\')  || strcmp(mnemonic_separator,'!') ...
            || strcmp(mnemonic_separator,'|') || strcmp(mnemonic_separator,'@')...
            || strcmp(mnemonic_separator, '/')
        mnemonic_startpos = findstr(char(fcsheader'),[mnemonic_name,mnemonic_separator]);
        if isempty(mnemonic_startpos)
            mneval = [];
            return;
        end
        mnemonic_length = length(mnemonic_name);
        mnemonic_stoppos = mnemonic_startpos + mnemonic_length;
        next_slashes = findstr(char(fcsheader(mnemonic_stoppos+1:end)'),mnemonic_separator);
        next_slash = next_slashes(1) + mnemonic_stoppos;

        mneval = char(fcsheader(mnemonic_stoppos+1:next_slash-1)');
    elseif strcmp(mnemonic_separator,'FF')
        mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
        if isempty(mnemonic_startpos)
            mneval = [];
            return;
        end
        mnemonic_length = length(mnemonic_name);
        mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
        next_formfeeds = find( fcsheader(mnemonic_stoppos+1:end) == 12);
        next_formfeed = next_formfeeds(1) + mnemonic_stoppos;

        mneval = char(fcsheader(mnemonic_stoppos + 1 : next_formfeed-1)');
    elseif strcmp(mnemonic_separator,'TAB')
        mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
        if isempty(mnemonic_startpos)
            mneval = [];
            return;
        end
        mnemonic_length = length(mnemonic_name);
        mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
        next_formfeeds = find( fcsheader(mnemonic_stoppos+1:end) == 9);
        next_formfeed = next_formfeeds(1) + mnemonic_stoppos;

        mneval = char(fcsheader(mnemonic_stoppos + 1 : next_formfeed-1)');

    elseif strcmp(mnemonic_separator, 'LF')
        mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
        if isempty(mnemonic_startpos)
            mneval = [];
            return;
        end
        mnemonic_length = length(mnemonic_name);
        mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
        next_linefeeds = find( fcsheader(mnemonic_stoppos+1:end) == 10);
        next_linefeed = next_linefeeds(1) + mnemonic_stoppos;

        mneval = char(fcsheader(mnemonic_stoppos + 1 : next_linefeed-1)');

    elseif strcmp(mnemonic_separator, 'RS')
        mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
        if isempty(mnemonic_startpos)
            mneval = [];
            return;
        end
        mnemonic_length = length(mnemonic_name);
        mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
        next_linefeeds = find( fcsheader(mnemonic_stoppos+1:end) == 30);
        next_linefeed = next_linefeeds(1) + mnemonic_stoppos;

        mneval = char(fcsheader(mnemonic_stoppos + 1 : next_linefeed-1)');
    end
end