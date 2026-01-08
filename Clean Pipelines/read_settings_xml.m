function[y_pos_table] = read_settings_xml(xml_file)
    
    doc = xmlread(xml_file);
    
    y_nodes = doc.getElementsByTagName('ELECTRODE_YPOS');
    
    assert(y_nodes.getLength > 0, 'No ELECTRODE_YPOS node found.');
    
    y_elem = y_nodes.item(0);
    
    % Read attributes like CH48="720"
    attributes = y_elem.getAttributes;
    n_channels = attributes.getLength;

    channel = nan(n_channels, 1);
    y_pos = nan(n_channels, 1);

    for i = 0:n_channels-1
        current_channel = attributes.item(i);
        channel_name = char(current_channel.getName);
        depth = str2double(current_channel.getValue);

        tok = regexp(channel_name, '^CH(\d+)$', 'tokens', 'once');

        if ~isempty(tok)
            channel(i+1) = str2double(tok{1});
            y_pos(i+1) = depth;
        end
    
    end
    y_pos_table = table(channel, y_pos, 'VariableNames',...
        {'ChannelID', 'y_pos'});

end