# -*- coding: utf-8 -*-

def load_read_data(fname, read_number=None, start=None, end=None, events=True, instance=0):
    """ Loads the read data for a specified read from a read fast5 file (taken from chimaera).
    :param fname: The fast5 file to read data from.
    :param read_number: Optional read number, for use with multi-read files.
    :param start: Start index of events to return. If None then starts at 0.
    :param end: End index of events to return. If None then returns to the end
        of the data.
    :param events: Flag indicating whether event data should be returned.
    :param instance: Specifies the 3 digit number of the EventDetection analyses group
        to pull events from.
    :returns: A tuple containing the data as a numpy record array, the channel,
        the read_number, and metadata from the fast5 file.
    """
    instance = str(instance).zfill(3)
    a = fname.rfind('_file')
    file_number = 0
    if a != -1:
        b = fname.rfind('_strand')
        if b != -1 and b > a:
            file_number = int(fname[a+5:b])
    with h5py.File(fname, 'r') as fast5:
        if 'file_version' not in fast5.attrs.keys():
            raise Exception('Cannot read from fast5 file. No version number specified.')
        if fast5.attrs['file_version'] != 1.0:
            raise Exception('Cannot read from fast5 file. Version not recognized.')
        if read_number is None:
            read_path = 'Analyses/EventDetection_{}/Reads'.format(instance)
            if read_path not in fast5:
                raise Exception('Could not parse fast5 file. No event data found.')
            reads_group = fast5[read_path]
            read_names = reads_group.keys()
            tokens = read_names[0].rsplit('_')
            read_number = int(tokens[1])
        channel_attrs = fast5['UniqueGlobalKey/channel_id'].attrs
        channel = int(channel_attrs['channel_number'])
        sampling_rate = float(channel_attrs['sampling_rate'])
        sample_size = 1.0 / sampling_rate
        read_path = 'Analyses/EventDetection_{}/Reads/Read_{}'.format(instance, read_number)
        if read_path not in fast5:
            raise Exception('Could not parse fast5 file. '
                            'No event data found for instance {} read {}.'.format(instance,
                                                                                  read_number))
        read_group = fast5[read_path]
        uuid = '00000000-0000-0000-0000-000000000000'
        if 'read_id' in read_group.attrs.keys():
            uuid = read_group.attrs['read_id']
        h5path = '{}/Events'.format(read_path)
        if events:
            event_data = fast5[h5path]
            if start is None:
                start = 0
            if end is None or end > event_data.size:
                end = event_data.size
            if event_data['start'].dtype.kind in ['i', 'u']:
                descr = [(x[0], 'float64') if x[0] in ('start', 'length') else x
                         for x in event_data.dtype.descr]
                with event_data.astype(np.dtype(descr)):
                    data = event_data[start:end]
                    data['start'] *= sample_size
                    data['length'] *= sample_size
            else:
                data = event_data[start:end]
        else:
            data = None
    meta = {'h5path': h5path, 'file_number': file_number, 'sampling_rate': sampling_rate,
            'uuid': uuid}
    return data, channel, read_number, meta
